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
While this API should be operational, it was not tested with any real BASIC assembly inputs.
Correctness of all assumptions is only guaranteed when the Link between part N and N+1 is the same for every construct.

In assembly_template.py:
- After all import statements, insert (if DNABOT and pipette_opt project folders are in the same folder):
    # Importing API from the neighbouring folder with pipette_opt project
    # Will be different in the future (maybe pipette_opt will be a package?)
    import sys
    sys.path.append('../pipette_opt')
    import pipette_opt as ppopt

- Lines 76-80
    Comment the whole 'part transfers' section
- Line 81:
    Insert: 
            action_list = ppopt.basic_part_transfer_actions(final_assembly_dict,1.5,pipette,method=...)
            ppopt.basic_execute(action_list,pipette,magbead_plate,destination_plate)
    cf. Start-Stop Assembly for choosing the method
"""

# -----------------------class definitions-----------------------------
# FOR OPTIMISATION
#operations
class Oper:
    def __init__(self, part, well):
        self.part = part
        self.well = well
        self.changed = False

    def __str__(self):  # for printing the subset's part type and wells out
        strRep = self.part + ' -> w' + str(self.well)
        return strRep

#operation subsets - needed for LP method
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
    # PART 1 Convert data into internal format

    # PART 1.1 Empty objects to fill
    w=[] # the array with information on constructs and parts, used as input to the algorithms
    dic = {'constructs': {}, 'parts': {}}  # the dictionary to decode algorithm outputs into action list
    address = {}  # addresses of part types in w
    reqvols ={} # list of required liquid volumes (needed for capacity)

    # PART 1.2 Declare auxiliary objects
    ignorelist = [] # list of parts to be ignored
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
    air_gap=1
    caps = capacities(reqvols=reqvols, pipcap=pipette.max_volume, airgap=air_gap)

    # PART 2 Solve the problem
    fin = [] # output operations (internal format)

    # PART 2.1 call algorithm specified by method
    if (method[0:2] == 'LP'): # LP-based
        if (len(method) == 2):
            tsp_method(w, fin, reord=None, filename=None, caps=caps)
        else:
            tsp_method(w, fin, method[3:], filename=None, caps=caps)
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
    cost = route_cost_with_w(fin, w, caps)
    savings = len(fin) - cost
    percentsavings = savings/len(fin)*100
    print('\npipette_opt: '+str(savings) + ' pipette tips saved (' + str(percentsavings) + '%)\n')

    # PART 2.3 Record in a .json file
    rec('Start-Stop',w,fin,dic,caps)

    # PART 3 Convert internal-output operations into action list

    # PART 3.1 Initialise action list, define unchanging parameters
    action_list = tuple()
    new_tip='once'


    # PART 3.2 Define auxiliary variables
    added = np.zeros((len(w), len(w[0])))  # tells which parts were added to which well
    cost = 1

    # PART 3.3 The first operation in fin
    added[fin[0].well][address[fin[0].part[0]]] = 1
    fin[0].changed = True
    part_source=dic['parts'][fin[0].part]['part_liqloc']
    part_vol = reqvols[fin[0].part]
    part_dest = dic['constructs'][fin[0].well]['con_liqloc']

    # PART 3.4 All later operations
    for i in range(1, len(fin)):
        # get operation cost
        cost = cost_func_with_w(fin[0:i], fin[i], w, added, caps)
        if (cost == 1):
            fin[i].changed = True
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


# ---------------------BASIC ASSEMBLY---------------------------
def basic_part_transfer_actions(method, final_assembly_dict, part_vol, pipette):
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
        action_list += basic_part_transfer_actions_onelen(method, one_dict, part_vol, pipette)

    return action_list


def basic_part_transfer_actions_onelen(method, final_assembly_dict, part_vol, pipette):
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
    caps = capacities(reqvols=reqvols, pipcap=pipette.max_volume, airgap=air_gap)

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
    cost = route_cost_with_w(fin, w, caps)
    savings = len(fin) - cost
    percentsavings = savings / len(fin) * 100
    print('\npipette_opt: ' + str(savings) + ' pipette tips saved (' + str(percentsavings) + '%)\n')

    # PART 2.3 record
    rec('BASIC', w, fin, dic, caps)

    # PART 3 Convert internal-output operations into an action list

    # PART 3.1 Initialise action list, define unchanging parameters
    action_list = tuple()
    new_tip = 'once'

    # PART 3.2 Define auxiliary variables
    added = np.zeros((len(w), len(w[0])))  # tells which parts were added to which well

    # PART 3.3 The first operation in fin
    added[fin[0].well][address[fin[0].part[0]]] = 1
    fin[0].changed = True
    part_source = dic['parts'][fin[0].part]['part_liqloc']
    part_vol = reqvols[fin[0].part]
    part_dest = [dic['constructs'][fin[0].well]['con_liqloc']]

    # PART 3.4 All later operations
    for i in range(1, len(fin)):
        # get operation cost
        cost = cost_func_with_w(fin[0:i], fin[i], w, added, caps)
        if (cost == 1):
            fin[i].changed = True
        added[fin[i].well][address[fin[i].part[0]]] = 1

        # act accroding to cost
        if (cost == 0):  # if 0, tip is unchanged
            part_dest = [dic['constructs'][fin[i].well]['con_liqloc']] + part_dest  # mind that we add at the beginning
        else:  # if 1, there is new tip, so...
            # record last tip's actions
            action_list += ((part_source, part_dest, part_vol, new_tip, air_gap),)

            # start recording new tip details
            part_source = dic['parts'][fin[i].part]['part_liqloc']
            part_vol = reqvols[fin[i].part]
            part_dest = [dic['constructs'][fin[i].well]['con_liqloc']]

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

        # PART 4 Drop used tip
        pipette.drop_tip()

    return


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