# PIPETTE_OPT.PY
# This is a 'Grand Central' program, which overlooks the algorithms.
# It acts as an adapter between assemblies and the algorithms, as well as providing a way to test the algorithms easily.
# Sort of like an API
# v0.0.1, 10.8.2020

from opentrons import instruments, labware

from auxil import *
from tsp_method import tsp_method
from statespace_methods import nns, greedy_tree

from vis_cont import rec

"""
THE START-STOP ASSEMBLY CAN NOW BENEFIT FROM THE PIPETTE TIP OPTIMISATION ALGORITHMS!
Just make the following changes to assembly.py...
(Here, the indents are used to aid legibility; please do appropriate formatting after inserting)

- After all import statements, insert (if dna_assembler and pipette_opt project folders are in the same folder):
    # Importing API from the neighbouring folder with pipette_opt project
    # Will be different in the future (maybe pipette_opt will be a package?)
    import sys
    sys.path.append('../pipette_opt')
    import pipette_opt as ppopt

- Line 323
    After:
        p300 = protocol.load_instrument('p300_single', 'right', tip_racks=[tiprack_1]) 
    insert:
        tiprack_2 = protocol.load_labware('opentrons_96_tiprack_10ul', 6)  # p10 pipette for kirill's tests
        p10 = protocol.load_instrument('p10_single', 'left', tip_racks=[tiprack_2])  # p10 pipette for kirill's tests

- Line 384
    Replace: total_mm = assembler.do_assembly(destination_location_manager=dest_loc_manager)
    by:      total_mm = ppopt.startstop_do_assembly(destination_location_manager=dest_loc_manager, dna_pipette=p10, method=...)
    
    ! method can be...
        ...'LP' (recommended, requires GUROBI optimiser)
        ...'LP+sametogether'
        ...'LP+leastout'
        ...'LP+greedy'
        ...'LP+nearest neighbour'
        ...'LP+iddfs depth 2'
        ...'Nearest Neighbour'
        ...'iddfs depth 2'
        ...'Greedy'
        ...'Nearest Neighbour+sametogether' (recommended Open Access option)
        ...'iddfs depth 2+sametogether'
        ...'Greedy+sametogether'
        Note that Hub-spoke is NOT supported
    
- Line 388
    Replace: assembler.distribute_dna(pipette=p300)
    by:      ppopt.startstop_distribute_dna(assembler,pipette=p10)
"""

"""
FOR BASIC ASSEMBLY

In dnabot_app.py:
- After all import statements, insert:
     # Importing API from the neighbouring folder with pipette_opt project
    # Will be different in the future (maybe pipette_opt will be a package?)
    import sys
    sys.path.append('../pipette_opt')
    import pipette_opt as ppopt
    
- Line 105 (after the 'calculate OT2 script variables' section). Insert (method selection the same as for Start-Stop):
    print('pipette_opt: optimising part distribution order...')
    ppopt_action_list = ppopt.basic_part_transfer_actions(method=..., final_assembly_dict=final_assembly_dict, part_vol=1.5, pipette_vol = 10)

- Line 117 - add a kwarg action_list=ppopt_action_list to get:
    generate_ot2_script(F_ASSEMBLY_FNAME, os.path.join(
        template_dir_path, F_ASSEMBLY_TEMP_FNAME),
        final_assembly_dict=final_assembly_dict,
        tiprack_num=final_assembly_tipracks,action_list=ppopt_action_list)
    


In assembly_template.py:
- In final_assembly() function:
    Comment the whole 'part transfers' section
    After it, put  a line: 
        basic_execute(action_list, pipette, magbead_plate, destination_plate)
        
- After the final_assembly() definition, insert the whole definition of basic_execute() from this file
"""

"""
FOR MOCLO ASSEMBLY

In moclo_transform_generator.py:
- After all import statements, insert:
    # Importing API from the neighbouring folder with pipette_opt project
    # Will be different in the future (maybe pipette_opt will be a package?)
    import sys
    sys.path.append('../pipette_opt')
    import pipette_opt as ppopt
    
- Line 45 (after calling the generate_and_save_output_plate_maps() function ). Insert (method selection the same as for Start-Stop):
    # PIPETTE_OPT
	# knowing that in this assembly a p10 pipette is used to distribute 2uL part aliquots, create action list
	ppopt_action_list = ppopt.moclo_part_transfer_actions(method=..., combinations_to_make=combinations_to_make, part_vol=2, pipette_vol=10)

- Line 117 - when calling the create_protocol() function, also include a ppopt_action_list=ppopt_action_list argument:
    create_protocol(dna_plate_map_dict, combinations_to_make, config['protocol_template_path'], config['output_folder_path'], ppopt_action_list=ppopt_action_list)

- Line 202 - when declaring the create_protocol() function, add an extra ppopt_action list argument:
    def create_protocol(dna_plate_map_dict, combinations_to_make, protocol_template_path, output_folder_path, ppopt_action_list):
    
- Line 212 (after writing dna_plate_dict and combinations_to_make in the script) - also write action_list in the script:
    # PIPETTE_OPT
		protocol_file.write('action_list = ' + json.dumps(ppopt_action_list) + '\n\n')
		

In moclo_transform_template.py
- Line 246 - COMMENT OUT the following section:
    combinations_by_part = {}
    for i in combinations_to_make:
        name = i["name"]
        for j in i["parts"]:
            if j in combinations_by_part.keys():
                combinations_by_part[j].append(name)
            else:
                combinations_by_part[j] = [name]
                
- Lines 255-278 (right after the just commented-out section):
    COMMENT OUT the whole section titled #This section of the code combines and mix the DNA parts according to the combination list
    
- Line 280 (right after that):
    insert the definition of the moclo_execute() function, which is given later in this file
    
- Line 311 (right after the moclo_execute() definition) - call moclo_execute():
    moclo_execute(action_list,p10_single)

"""


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
    letcodes = 'prctbadefghijklmnoqsuvxyz'  # one-letter ids standing in for part types
    i_letcodes = 0  # points to one-letter id of the (potential) next part type to record
    i_addr = 0 # points to address in w of the (potential) next part type to record

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
                parttype[part] = letcodes[i_letcodes]
                partnum[letcodes[i_letcodes]] = 0
                address[letcodes[i_letcodes]] = i_addr
                i_letcodes += 1
                i_addr += 1

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
                newcode = parttype[part] + str(partnum[parttype[part]])  # determine which entry to record
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
            tsp_method(w, fin, reord=None, caps=caps)
        else:
            tsp_method(w, fin, method[3:], caps=caps)
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


# ---------------------BASIC ASSEMBLY---------------------------
def basic_part_transfer_actions(method, final_assembly_dict, part_vol, pipette_vol):
    # PART 1: Group all the constructs by lengths
    assembly_dict_by_lens={}
    for key in final_assembly_dict.keys():
        # check if such length is already present
        match = False
        final_len = len(final_assembly_dict[key])
        for key_by_lens in assembly_dict_by_lens.keys():
            if(final_len == key_by_lens):
                match = True
                break

        # if yes, add
        if(match):
            assembly_dict_by_lens[key_by_lens][key] = final_assembly_dict[key]
        # if no, create a new entry with this length and this construct standing for it
        else:
            assembly_dict_by_lens[final_len] = {key: final_assembly_dict[key]}

    # PART 2: Make the action list
    # initialise
    action_list = tuple()
    # add actions for each group of subsets
    for one_dict in assembly_dict_by_lens.values():
        action_list += basic_part_transfer_actions_onelen(method, one_dict, part_vol, pipette_vol)

    return action_list


def basic_part_transfer_actions_onelen(method, final_assembly_dict, part_vol, pipette_volume):
    # PART 1 Convert data into internal format

    # PART 1.1 Empty objects to fill
    w = []  # the array with information on constructs and parts, used as input to the algorithms
    dic = {'constructs': {}, 'parts': {}}  # the dictionary to decode algorithm outputs into action list
    address = {}  # addresses of part types in w
    reqvols = {}  # list of required liquid volumes (needed for capacity)

    # PART 1.2 Declare auxiliary objects
    ignorelist =[] # !! what would we ignore with BASIC assembly?
    parttype = {}  # indices of parts to be recorded in the subsets list
    partnum = {}  # current number of each part type different species
    wellno = 0  # well counter
    letcodes = 'abcdefghijklmnopqrstuvxyz'  # one-letter ids standing in for part types
    i_letcodes = 0  # points to one-letter id of the (potential) next part type to record
    i_addr = 0  # points to address in w of the (potential) next part type to record


    # PART 1.3 Read the constructs
    for construct in final_assembly_dict.keys():
        # match well number in w with its location on the plate and construct
        dic['constructs'][wellno] = {'con_liqloc': construct}

        # get parts
        parts = final_assembly_dict[construct] # get list of parts of the cinstruct
        well_parts=[] # will store all codes of parts of current well
        for i in range(0,len(parts)):
            # skip an ignored part type
            ig = False
            for ignore in ignorelist:
                if (parts[i] == ignore):
                    ig = True
                    break
            if(ig):
                continue

            #  if this is the first entry, fill the address dictionary
            if (wellno == 0):
                parttype[i] = letcodes[i_letcodes]
                partnum[letcodes[i_letcodes]] = 0
                address[letcodes[i_letcodes]] = i_addr
                i_letcodes += 1
                i_addr += 1

            # determine part name
            partloc = parts[i]

            # check if such location is already in the dictionary
            match = False
            for partcode in dic['parts'].keys():
                if (partloc == dic['parts'][partcode]['part_liqloc']):  # yes, record the part code in w
                    well_parts.append(partcode)
                    match = True
                    break

            if (not match):  # if no, update the dictionary
                newcode = parttype[i] + str(partnum[parttype[i]])  # determine which entry to record
                dic['parts'][newcode] = {'part_liqloc': parts[i]}  # put the entry into dictionary
                well_parts.append(newcode)  # record part code in w
                reqvols[newcode] = part_vol  # record the required volume for the part - which is the same in BASIC
                partnum[parttype[i]] += 1  # update number of parts of this class

        w.append(well_parts.copy())  # record parts of the latest well
        wellno += 1  # proceeding to next well

    # PART 1.4 Get capacities
    air_gap=1
    #if(type(pipette)!=)
    caps = capacities(reqvols=reqvols, pipcap=pipette_volume, airgap=air_gap)

    # PART 2 Solve the problem
    fin = []  # output operations (internal format)

    # PART 2.1 call algorithm specified by method
    if (method[0:2] == 'LP'):  # LP-based
        if (len(method) == 2):
            tsp_method(w, fin, reord=None, caps=caps)
        else:
            tsp_method(w, fin, method[3:], caps=caps)
    else:  # state-space search
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
    percentsavings = savings / len(fin) * 100
    print('pipette_opt:\n')
    print(str(savings) + ' pipette tips saved (' + str(percentsavings) + '%)')
    print('Thus ' + str(cost) + ' tips required\n')

    # PART 2.3 record
    rec('BASIC', w, fin, dic, caps)

    # PART 3 Convert internal-output operations into an action list

    # PART 3.1 Initialise action list, define unchanging parameters
    action_list = tuple()
    new_tip = 'once'

    # PART 3.2 The first operation in fin
    part_source = dic['parts'][fin[0].part]['part_liqloc']
    part_vol = reqvols[fin[0].part]
    part_dest = [dic['constructs'][fin[0].well]['con_liqloc']].copy()

    # PART 3.4 All later operations
    for i in range(1, len(fin)):
        # act accroding to cost
        if not (fin[i].changed):  # if the tip is unchanged
            part_dest += [dic['constructs'][fin[i].well]['con_liqloc']]
        else:  # if there is a new tip...
            # record last tip's actions
            action_list += ((part_source, part_dest, part_vol, new_tip, air_gap),)

            # start recording new tip details
            part_source = dic['parts'][fin[i].part]['part_liqloc']
            part_vol = reqvols[fin[i].part]
            part_dest = [dic['constructs'][fin[i].well]['con_liqloc']].copy()

        # if this is the end, just record last tip's actions
        if (i == len(fin) - 1):
            action_list += ((part_source, part_dest, part_vol, new_tip, air_gap),)

    return action_list


def basic_execute(action_list, pipette, magbead_plate, destination_plate):
    for action in action_list:
        # PART 1 Get directions
        source_well, destination_wells, volume, new_tip, air_gap = action

        # PART 2 Get new tip
        pipette.pick_up_tip()

        # PART 3 Aspirate
        # 3a) with air gaps
        if(air_gap != 0):
            for i in range(0, len(destination_wells) - 1):
                pipette.aspirate(volume, magbead_plate.wells(source_well))
                pipette.air_gap(air_gap)
            pipette.aspirate(volume, magbead_plate.wells(source_well))
        # 3b) no air gaps
        else:
            pipette.aspirate(volume*len(destination_wells),magbead_plate.wells(source_well))


        # PART 4 Distribute
        for i in range(0, len(destination_wells) - 1):
            pipette.dispense(volume + air_gap, destination_plate.wells(destination_wells[i]))
        pipette.dispense(volume, destination_plate.wells(destination_wells[- 1]))

        # PART 5 Drop used tip
        pipette.drop_tip()

    return


# ---------------------MOCLO ASSEMBLY---------------------------
def moclo_part_transfer_actions(method, combinations_to_make, part_vol, pipette_vol):
    # PART 1: Group all the constructs by lengths
    assembly_dict_by_lens={}
    for comb in combinations_to_make:
        # check if such length is already present
        match = False
        final_len = len(comb['parts'])
        for key_by_lens in assembly_dict_by_lens.keys():
            if(final_len == key_by_lens):
                match = True
                break

        # if yes, add
        if(match):
            assembly_dict_by_lens[key_by_lens].append(comb)
        # if no, create a new entry with this length and this construct standing for it
        else:
            assembly_dict_by_lens[final_len] = [comb]

    # PART 2: Make the action list
    # initialise
    action_list = tuple()
    # add actions for each group of subsets
    for one_dict in assembly_dict_by_lens.values():
        action_list += moclo_part_transfer_actions_onelen(method, one_dict, part_vol, pipette_vol)

    return action_list


def moclo_part_transfer_actions_onelen(method, constructs, part_vol, pipette_volume):
    # PART 1 Convert data into internal format

    # PART 1.1 Empty objects to fill
    w = []  # the array with information on constructs and parts, used as input to the algorithms
    dic = {'constructs': {}, 'parts': {}}  # the dictionary to decode algorithm outputs into action list
    address = {}  # addresses of part types in w
    reqvols = {}  # list of required liquid volumes (needed for capacity)

    # PART 1.2 Declare auxiliary objects
    ignorelist =[] # !! what would we ignore with BASIC assembly?
    parttype = {}  # indices of parts to be recorded in the subsets list
    partnum = {}  # current number of each part type different species
    wellno = 0  # well counter
    letcodes = 'abcdefghijklmnopqrstuvxyz'  # one-letter ids standing in for part types
    i_letcodes = 0  # points to one-letter id of the (potential) next part type to record
    i_addr = 0  # points to address in w of the (potential) next part type to record


    # PART 1.3 Read the constructs
    for construct in constructs:
        # match well number in w with its location on the plate and construct
        dic['constructs'][wellno] = {'con_name': construct['name']}

        # get parts
        parts = construct['parts'] # get list of parts of the cinstruct
        well_parts=[] # will store all codes of parts of current well
        for i in range(0,len(parts)):
            # skip an ignored part type
            ig = False
            for ignore in ignorelist:
                if (parts[i] == ignore):
                    ig = True
                    break
            if(ig):
                continue

            #  if this is the first entry, fill the address dictionary
            if (wellno == 0):
                parttype[i] = letcodes[i_letcodes]
                partnum[letcodes[i_letcodes]] = 0
                address[letcodes[i_letcodes]] = i_addr
                i_letcodes += 1
                i_addr += 1

            # determine part name
            part = parts[i]

            # check if such location is already in the dictionary
            match = False
            for partcode in dic['parts'].keys():
                if (part == dic['parts'][partcode]['part_name']):  # yes, record the part code in w
                    well_parts.append(partcode)
                    match = True
                    break

            if (not match):  # if no, update the dictionary
                newcode = parttype[i] + str(partnum[parttype[i]])  # determine which entry to record
                dic['parts'][newcode] = {'part_name': parts[i]}  # put the entry into dictionary
                well_parts.append(newcode)  # record part code in w
                reqvols[newcode] = part_vol  # record the required volume for the part - which is the same in BASIC
                partnum[parttype[i]] += 1  # update number of parts of this class

        w.append(well_parts.copy())  # record parts of the latest well
        wellno += 1  # proceeding to next well

    # PART 1.4 Get capacities
    air_gap=1
    #if(type(pipette)!=)
    caps = capacities(reqvols=reqvols, pipcap=pipette_volume, airgap=air_gap)

    # PART 2 Solve the problem
    fin = []  # output operations (internal format)

    # PART 2.1 call algorithm specified by method
    if (method[0:2] == 'LP'):  # LP-based
        if (len(method) == 2):
            tsp_method(w, fin, reord=None, caps=caps)
        else:
            tsp_method(w, fin, method[3:], caps=caps)
    else:  # state-space search
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
    percentsavings = savings / len(fin) * 100
    print('pipette_opt:\n')
    print(str(savings) + ' pipette tips saved (' + str(percentsavings) + '%)')
    print('Thus ' + str(cost) + ' tips required\n')

    # PART 3 Convert internal-output operations into an action list

    # PART 3.1 Initialise action list, define unchanging parameters
    action_list = tuple()
    new_tip = 'once'

    # PART 3.2 The first operation in fin
    part_source = dic['parts'][fin[0].part]['part_name']
    part_vol = reqvols[fin[0].part]
    part_dest = [dic['constructs'][fin[0].well]['con_name']].copy()

    # PART 3.3 All later operations
    for i in range(1, len(fin)):
        # act according to cost
        if not (fin[i].changed):  # if 0, tip is unchanged
            part_dest += [dic['constructs'][fin[i].well]['con_name']]
        else:  # if there is a new tip...
            # record last tip's actions
            action_list += ((part_source, part_dest, part_vol, new_tip, air_gap),)

            # start recording new tip details
            part_source = dic['parts'][fin[i].part]['part_name']
            part_vol = reqvols[fin[i].part]
            part_dest = [dic['constructs'][fin[i].well]['con_name']].copy()

        # if this is the end, just record last tip's actions
        if (i == len(fin) - 1):
            action_list += ((part_source, part_dest, part_vol, new_tip, air_gap),)

    return action_list

'''
def moclo_execute(action_list, pipette):
    for action in action_list:
        # PART 1 Get directions
        source_well, destination_wells, volume, new_tip, air_gap = action

        # PART 2 Get new tip
        pipette.pick_up_tip()

        # PART 3 Aspirate
        # 3a) with air gaps
        if(air_gap != 0):
            for i in range(0, len(destination_wells) - 1):
                pipette.aspirate(volume, find_dna(source_well,dna_plate_map_dict,dna_plate_dict))
                pipette.air_gap(air_gap)
            pipette.aspirate(volume, find_dna(source_well,dna_plate_map_dict,dna_plate_dict))
        # 3b) no air gaps
        else:
            pipette.aspirate(volume*len(destination_wells),find_dna(source_well,dna_plate_map_dict,dna_plate_dict))


        # PART 4 Distribute
        for i in range(0, len(destination_wells) - 1):
            pipette.dispense(volume + air_gap, find_combination(destination_wells[i],combinations_to_make))
        pipette.dispense(volume, find_combination(destination_wells[-1],combinations_to_make))

        # PART 5 Drop used tip
        pipette.drop_tip()

    return
'''

# -------------------------- MAIN (TEST ONLY) -------------------------
def main():
    final_assembly_dict = {"A1": ["A7", "B7", "C7", "D7", "E7"], "B1": ["F7", "B7", "G7", "H7", "E7"],
                           "C1": ["A7", "A8", "B8", "C8", "E7"], "D1": ["F7", "A8", "G7", "D7", "E7"]}

    PIPETTE_MOUNT = 'right'
    MAG_PLATE_TYPE = 'corning_96_wellplate_360ul_flat'
    MAG_PLATE_POSITION = '1'
    DESTINATION_PLATE_TYPE = 'opentrons_96_tiprack_300ul'
    TEMPDECK_SLOT = '4'

    pipette = instruments.P10_Single(mount=PIPETTE_MOUNT)
    magbead_plate = labware.load(MAG_PLATE_TYPE, MAG_PLATE_POSITION)
    destination_plate = labware.load(DESTINATION_PLATE_TYPE, TEMPDECK_SLOT, share=True)

    action_list = basic_part_transfer_actions('LP+sametogether',final_assembly_dict,1.5,pipette)
    #basic_execute(action_list,pipette,magbead_plate,destination_plate)

if __name__ == "__main__":
    main()