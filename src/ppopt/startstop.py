# startstop - API for the 'github.com/zoltuz/dna_assembler' automated assmembly package (Start-Stop standard)
# v0.1.0, 31.5.21

"""
Make the following changes to assembly.py:
(Here, the indents are used to aid legibility; appropriate formatting after inserting is required)

- After all import statements, insert (if dna_assembler and pipette_opt project folders are in the same folder):
    # Importing API for ppopt optimisation of pipette tip consumption
    import pipette_opt as ppopt

- Line 334
    After:
        p300 = protocol.load_instrument('p300_single', 'right', tip_racks=[tiprack_1])
    insert:
        tiprack_2 = protocol.load_labware('opentrons_96_tiprack_10ul', 6)  # for example, use p10 pipette
        p10 = protocol.load_instrument('p10_single', 'left', tip_racks=[tiprack_2])  # for example, use p10 pipette

- Line 484
    Replace: total_mm = assembler.do_assembly(destination_location_manager=dest_loc_manager)
    by:      total_mm = startstop_do_assembly(assembly,destination_location_manager=dest_loc_manager, dna_pipette=p10, method=...)

    ! method can be...
        ...'LP' (recommended, requires GUROBI optimiser)
        ...'LP+sametogether'
        ...'LP+leastout'
        ...'LP+greedy'
        ...'LP+nearest neighbour'
        ...'LP+iddfs depth 2'
        ...'Nearest Neighbour'
        ...'nns depth 2'
        ...'Greedy'
        ...'Nearest Neighbour+sametogether' (recommended Open Access option)
        ...'nns depth 2+sametogether'
        ...'Greedy+sametogether'
        Note that Hub-spoke is NOT supported

- Line 487
    Replace: assembler.distribute_dna(pipette=p300)
    by:      startstop_distribute_dna(assembler,pipette=p10)
"""

from ppopt.statespace import nns, greedy_tree
from ppopt.lp import lp_method
from ppopt.dp import dp_method
from ppopt.auxil import *

# ---------------------START-STOP ASSEMBLY---------------------
# Modified LevelZeroAssembler.do_assembly() function
# Also needs a DNA-distributing pipette and opt method as extra inputs
# Now producing an OPTIMISED action list
def startstop_do_assembly(assembly, destination_location_manager, dna_pipette, method):
    # reset action list
    assembly.action_list = list()
    assembly.ddh2o_list = list()
    total_ddH2O_volume = 0.0
    for construct in assembly.current_constructs:
        # get a location for the levelzero construct from the location manager
        dest_location = destination_location_manager.get_a_location()
        construct.liquid_location = dest_location
        total_ddH2O_volume += construct.reaction_components['ddH2O']
        assembly.ddh2o_list.append((dest_location, construct.reaction_components['ddH2O']))
    assembly.action_list = startstop_actionlist(assembly, method, dna_pipette)
    return total_ddH2O_volume


# Creates an optimised action list
def startstop_actionlist(assembly,method, pipette):
    # PART 0: handle the vector backbones
    backbone_action_list, ignorelist = startstop_handle_backbones(assembly,method, pipette)


    # PART 1 Convert data into internal format

    # PART 1.1 Empty objects to fill
    w=[] # the array with information on constructs and parts, used as input to the algorithms
    dic = {'constructs': {}, 'parts': {}}  # the dictionary to decode algorithm outputs into action list
    address = {}  # addresses of part types in w
    reqvols ={} # list of required liquid volumes (needed for capacity)

    # PART 1.2 Declare auxiliary objects
    parttype = {}  # indices of parts to be recorded in the subsets list
    partnum = {}  # current number of each part type different species
    wellno = 0  # well counter
    partpos=0 # position of a given part type in the array w

    # PART 1.3 Read the constructs
    for construct in assembly.current_constructs:
        # match well number in w with its location on the plate and construct
        dic['constructs'][wellno] = {'con_name': construct.construct_name,
                                     'con_liqloc': construct.liquid_location}

        # get parts
        parts = construct.plasmid # get list of parts of the construct
        well_parts=[] # will store all codes of parts of current well
        for part in parts.keys():
            # skip an ignored part type
            ig = False
            for ignore in ignorelist:
                if (part == ignore):
                    ig = True
                    break
            if(ig):
                continue

            #  if this is the first entry, fill the address dictionary
            if (wellno == 0):
                parttype[part] = partpos
                partnum[partpos] = 0
                partpos+=1

            # determine part name
            partname = parts[part].name

            # check if such name is already in the dictionary
            match = False
            for partcode in dic['parts'].keys():
                if (partname == dic['parts'][partcode]['part_name']):  # yes, record the part code in w
                    well_parts.append(partcode)
                    match = True
                    break

            if (not match):  # if no, update the dictionary
                newcode = (parttype[part],partnum[parttype[part]])  # determine which entry to record
                dic['parts'][newcode] = {'part_name': partname,
                                         'part_liqloc': parts[part].liquid_location}  # put the entry into dictionary
                well_parts.append(newcode)  # record part code in w
                reqvols[newcode] = parts[part].required_volume  # record the required volume for the part
                partnum[parttype[part]] += 1  # update number of parts of this class

        w.append(well_parts.copy())  # record parts of the latest well
        wellno += 1  # proceeding to next well

    # PART 1.4 Get capacities
    air_gap = 1
    caps = capacities(reqvols=reqvols, pipcap=pipette.max_volume, airgap=air_gap)

    # PART 2 Solve the problem
    fin = [] # output operations (internal format)

    # PART 2.1 call algorithm specified by method
    if (method[0:2] == 'LP'): # LP-based
        if (len(method) == 2):
            lp_method(w, fin, reord=None, caps=caps, maxtime=1)
        else:
            lp_method(w, fin, method[3:], caps=caps, maxtime=1)
    elif (method[0:2] == 'DP'):  # LP-based
        if (len(method) == 2):
            dp_method(w, fin, reord=None, caps=caps)
        else:
            dp_method(w, fin, reord=method[3:], caps=caps)
    else: # state-space search
        # determine reordeing
        if (method[-12:] == 'sametogether'):
            reord = 'sametogether'
        else:
            reord = None
        # get solution
        if (method[:7] == 'Nearest'):
            nns(w, fin, 1, reord, caps)
        elif (method[:3] == 'nns'):
            nns(w, fin, 2, reord, caps)
        elif (method[:6] == 'Greedy'):
            greedy_tree(w, fin, 'optimistic+cap', reord, caps)

    # PART 2.2 Print report on solution benefits
    cost = route_cost(fin)
    savings = len(fin) - cost
    percentsavings = savings/len(fin)*100
    print('\npipette_opt: '+str(savings) + ' pipette tips saved (' + str(percentsavings) + '%)')
    print('Thus ' + str(cost) + ' tips required\n')

    # PART 2.3 Record in a .p file
    rec('Start-Stop',w,fin,dic,caps)

    # PART 3 Convert internal-output operations into action list

    # PART 3.1 Initialise action list with the output of backbone handler, define unchanging parameters
    action_list = backbone_action_list
    new_tip = 'once'

    # PART 3.2 The first operation in fin
    part_source=dic['parts'][fin[0].part]['part_liqloc']
    part_vol = reqvols[fin[0].part]
    part_dest = dic['constructs'][fin[0].well]['con_liqloc'].copy()

    # PART 3.3 All later operations
    for i in range(1, len(fin)):
        # act accroding to whtehre the tip was changed

        if not (fin[i].changed): # if the tip is unchanged
            part_dest += dic['constructs'][fin[i].well]['con_liqloc']
        else: # if there is a new tip...
            # record last tip's actions
            action_list+=((part_source,part_dest,part_vol,new_tip,air_gap),)

            # start recording new tip details
            part_source = dic['parts'][fin[i].part]['part_liqloc']
            part_vol = reqvols[fin[i].part]
            part_dest = dic['constructs'][fin[i].well]['con_liqloc'].copy()

        # if this is the end, just record last tip's actions
        if(i == len(fin)-1):
            action_list += ((part_source, part_dest, part_vol, new_tip, air_gap),)

    return action_list


def startstop_distribute_dna(assembly, pipette):
    for action in assembly.action_list:
        # PART 1 Get directions
        source_well, destination_wells, volume, new_tip, air_gap = action

        # PART 2 Get new tip
        pipette.pick_up_tip()

        # PART 3 Aspirate
        # 3a) with air gaps
        if (air_gap != 0):
            for i in range(0, len(destination_wells) - 1):
                pipette.aspirate(volume, source_well)
                pipette.air_gap(air_gap)
            pipette.aspirate(volume, source_well)
        # 3b) no air gaps
        else:
            pipette.aspirate(volume * len(destination_wells), source_well)

        # PART 4 Distribute
        for i in range(0,len(destination_wells)-1):
            pipette.dispense(volume+air_gap, destination_wells[i])
        pipette.dispense(volume, destination_wells[len(destination_wells)-1])

        # PART 4 Drop used tip
        pipette.drop_tip()

    return


# handling the backbone vectors
def startstop_handle_backbones(assembly,method, pipette):
    """
    - In the overwhelming majority of cases, the vector backbones are all the same, so it's most efficient not
    to handle them as a DNA part but just distribute them prior to all parts.
    The returned backbone_action_list will go before all other actions, so the backbones are distributed before
    everything else. The returned ignorelist will make the algorithms ignore the backbones when reading the inputs

    - In the unlikely event that there are different backbones, they can just be condidered DNA parts. The returned
    EMPTY backbone_action_list and ignorelist will ensure that.
    """
    # PART 1: initial preparations

    # PART 1.1: initialise output lists
    backbone_action_list = tuple() # pipette actions to distribute tha backbones
    ignorelist = [] # will make the algorithms ignore backbones later, if necessary

    # PART 1.1: get the list of constructs
    constructs = assembly.current_constructs


    # PART 2: determine if all backbones are the same
    allsame = True
    for i in range(1, len(constructs)):
        if (constructs[i].plasmid['vector_backbone'].name != constructs[0].plasmid['vector_backbone'].name):
            allsame = False
            break

    # PART 3: if yes, fill action and ignore lists
    if (allsame):
        # PART 3.1: update the ignore list
        ignorelist.append('vector_backbone')

        # PART 3.2: fill the action list
        # PART 3.2.1: get variables common for all actions
        bb = constructs[0].plasmid['vector_backbone']
        bb_cap = capac(pipcap=pipette.max_volume, dose=bb.required_volume, airgap=0)
        part_source = bb.liquid_location
        part_vol = bb.required_volume
        new_tip = 'once'
        air_gap = 0

        # PART 3.2.2: get destination of the first construct
        part_dest = constructs[0].liquid_location.copy()

        # PART 3.2.3: get all other destinations
        for i in range(1, len(constructs)):
            # if no need to change the tip, just add the well to destinations
            if (((i % bb_cap) != 0)):
                part_dest += constructs[i].liquid_location
            # otherwise, put previous destinations into action list and start the new destination list
            else:
                backbone_action_list += ((part_source, part_dest, part_vol, new_tip, air_gap),)
                part_dest = constructs[i].liquid_location

        # record the last actions that remain unrecorded
        backbone_action_list += ((part_source, part_dest, part_vol, new_tip, air_gap),)


    # PART 4: return
    return backbone_action_list, ignorelist
