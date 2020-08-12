# PIPETTE_OPT.PY
# This is a 'Grand Central' program, which overlooks the algorithms.
# It acts as an adapter between assemblies and the algorithms, as well as providing a way to test the algorithms easily.
# Sort of like an API
# v0.0.1, 10.8.2020

from auxil import *
from tsp_method import tsp_method
from statespace_methods import iddfs, greedy_tree

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
        ...'TSP'
        ...'TSP+sametogether'
        ...'TSP+leastout'
        ...'TSP+greedy'
        ...'TSP+nearest neighbour'
        ...'TSP+iddfs depth 2'
        ...'Nearest Neighbour'
        ...'iddfs depth 2'
        ...'Greedy'
        ...'Nearest Neighbour+sametogether'
        ...'iddfs depth 2+sametogether'
        ...'Greedy+sametogether'
        Note that Hub-spoke is NOT supported
    
- Line 388
    Replace: assembler.distribute_dna(pipette=p300)
    by:      ppopt.startstop_distribute_dna(assembler,pipette=p10)
"""

# -----------------------class definitions-----------------------------
#operations
class Oper:
    def __init__(self, part, well):
        self.part = part
        self.well = well

    def __str__(self):  # for printing the subset's part type and wells out
        strRep = self.part + ' -> w' + str(self.well)
        return strRep

#operation subsets - needed for TSP method
class Ss:
    def __init__(self, part, wellno):  # initialisation
        self.part = part
        self.wells = [wellno]

    def nuwell(self, wellno):  # record new well in the subset
        self.wells.append(wellno)

    def __str__(self):  # for printing the subset's part type and wells out
        strRep = self.part + '|'
        for i in range(0, len(self.wells)):
            strRep = strRep + ' ' + str(self.wells[i])
        return strRep

# ---------------------START-STOP ASSEMBLY---------------------
# Modified LevelZeroAssembler.do_assembly() function
# Also needs a DNA-distributing pipette and opt method as extra inputs
# Now producing an OPTIMISED action list
def startstop_do_assembly(self, destination_location_manager: LocationManager, dna_pipette, method):
    # reset action list
    self.action_list = list()
    self.ddh2o_list = list()
    total_ddH2O_volume = 0.0
    for construct in self.current_constructs:
        # get a location for the levelzero construct from the location manager
        dest_location = destination_location_manager.get_a_location()
        construct.liquid_location = dest_location
        total_ddH2O_volume += construct.reaction_components['ddH2O']
        self.ddh2o_list.append((dest_location, construct.reaction_components['ddH2O']))
    self.action_list = ppopt.startstop_actionlist(self, method, dna_pipette)
    return total_ddH2O_volume


# Creates an optimised action list
def startstop_actionlist(assembly,method, pipette):
    # PART 1 Convert data into internal format

    # PART 1.1 Empty objects to fill
    w=[] # the array with information on constructs and parts, used as input to the algorithms
    dic = {'constructs': {}, 'parts': {}}  # the dictionary to decode algorithm outputs into action list
    address = {}  # addresses of part types in w
    reqvols ={} # list of required liquid volumes (needed for capacity)

    # PART 1.2 Declare auxiliary objects
    ignorelist = ['vector_backbone'] # backbone part must be ingored for start-stop assembly
    parttype = {}  # indices of parts to be recorded in the subsets list
    partnum = {}  # current number of each part type different species
    wellno = 0  # well counter
    letcodes = 'prctabdefghijklmnoqsuvwxyz'  # one-letter ids standing in for part types
    i_letcodes = 0  # points to one-letter id of the (potential) next part type to record
    i_addr = 0 # points to address in w of the (potential) next part type to record

    # PART 1.3 Read the constructs
    for construct in assembly.current_constructs:
        # match well number in w with its location on the plate and construct
        dic['constructs'][wellno] = {'con_name': construct.construct_name,
                                     'con_liqloc': construct.liquid_location}

        # get parts
        parts = construct.plasmid # get list of parts of the cinstruct
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
    air_gap=1
    caps = capacities(reqvols=reqvols, pipcap=pipette.max_volume, airgap=air_gap)

    # PART 2 Solve the problem
    fin = [] # output operations (internal format)

    # PART 2.1 call algorithm specified by method
    if (method[0:3] == 'TSP'): # TSP-based
        if (len(method) == 3):
            tsp_method(w, fin, reord=None, filename=None, caps=caps)
        else:
            tsp_method(w, fin, method[4:], filename=None, caps=caps)
    else: # state-space search
        # determine reordeing
        if (method[-12:] == 'sametogether'):
            reord = 'sametogether'
        else:
            reord = None
        # get solution
        if (method[:7] == 'Nearest'):
            iddfs(w, fin, 1, True, reord, caps)
        elif (method[:5] == 'iddfs'):
            iddfs(w, fin, 2, True, reord, caps)
        elif (method[:6] == 'Greedy'):
            greedy_tree(w, fin, 'optimistic+cap', reord, caps)

    # PART 2.2 Print report on solution benefits
    cost = route_cost_with_w(fin, w, caps)
    savings = len(fin) - cost
    percentsavings = savings/len(fin)*100
    print('\npipette_opt: '+str(savings) + ' pipette tips saved (' + str(percentsavings) + '%)\n')

    # PART 3 Convert internal-output operations into action list

    # PART 3.1 Initialise action list, define unchanging parameters
    action_list = tuple()
    new_tip='once'


    # PART 3.2 Define auxiliary variables
    added = np.zeros((len(w), len(w[0])))  # tells which parts were added to which well
    cost = 1

    # PART 3.3 The first operation in fin
    added[fin[0].well][address[fin[0].part[0]]] = 1
    part_source=dic['parts'][fin[0].part]['part_liqloc']
    part_vol = reqvols[fin[0].part]
    part_dest = dic['constructs'][fin[0].well]['con_liqloc']

    # PART 3.4 All later operations
    for i in range(1, len(fin)):
        # get operation cost
        cost = cost_func_with_w(fin[0:i], fin[i], w, added, caps)
        added[fin[i].well][address[fin[i].part[0]]] = 1

        # act accroding to cost
        if(cost==0): # if 0, tip is unchanged
            part_dest = dic['constructs'][fin[i].well]['con_liqloc'] + part_dest # mind that we add at the beginning
        else: # if 1, there is new tip, so...
            # record last tip's actions
            action_list+=((part_source,part_dest,part_vol,new_tip,air_gap),)

            # start recording new tip details
            part_source = dic['parts'][fin[i].part]['part_liqloc']
            part_vol = reqvols[fin[i].part]
            part_dest = dic['constructs'][fin[i].well]['con_liqloc']

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
        for i in range(0,len(destination_wells)-1):
            pipette.aspirate(volume, source_well).air_gap(air_gap)
        pipette.aspirate(volume,source_well)

        # PART 4 Distribute
        for i in range(0,len(destination_wells)-1):
            pipette.dispense(volume+air_gap, destination_wells[i])
        pipette.dispense(volume, destination_wells[i])

        # PART 4 Drop used tip
        pipette.drop_tip()

    return
