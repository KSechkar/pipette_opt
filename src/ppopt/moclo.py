# moclo - API for the DNA-BOT automated assmembly package (Start-Stop standard)
# v0.1.0, 31.5.21

"""
FOR MOCLO ASSEMBLY

In moclo_transform_generator.py:
- After all import statements, insert:
    # Importing API for ppopt optimisation of pipette tip consumption
    from ppopt.moclo import moclo_part_transfer_actions

- Line 40 (after calling the generate_and_save_output_plate_maps() function ). Insert (method selection the same as for Start-Stop):
    # PPOPT: knowing that in this assembly a p10 pipette is used to distribute 2uL part aliquots, create action list
	ppopt_action_list = moclo_part_transfer_actions(method='LP', combinations_to_make=combinations_to_make, part_vol=2, pipette_vol=10)

- Line 44 - when calling the create_protocol() function, also include a ppopt_action_list=ppopt_action_list argument:
    create_protocol(dna_plate_map_dict, combinations_to_make, config['protocol_template_path'], config['output_folder_path'], ppopt_action_list=ppopt_action_list)

- Line 197 - when declaring the create_protocol() function, add an extra ppopt_action_list argument:
    def create_protocol(dna_plate_map_dict, combinations_to_make, protocol_template_path, output_folder_path, ppopt_action_list):

- Line 206 (after writing dna_plate_dict and combinations_to_make in the script) - also write action_list in the script:
    # PPOPT: record optimised action list
	protocol_file.write('action_list = ' + json.dumps(ppopt_action_list) + '\n\n')

In moclo_transform_template.py
- After all import statements, insert:
    # Importing API for ppopt optimisation of pipette tip consumption
    from ppopt.moclo import moclo_execute


- Line 268 - COMMENT OUT the following section:
    combinations_by_part = {}
    for i in combinations_to_make:
        name = i["name"]
        for j in i["parts"]:
            if j in combinations_by_part.keys():
                combinations_by_part[j].append(name)
            else:
                combinations_by_part[j] = [name]
- Lines 275-298 (right after the just commented-out section):
    COMMENT OUT the whole section titled #This section of the code combines and mix the DNA parts according to the combination list

- Line 300 (right after the just commented-out section) - insert definition of moclo_execute

- Line 330 (right after newly inserted defintion) - insert:
    # PPOPT: execute optimised list of pipette actions
    moclo_execute(action_list,p10_single)
"""

from ppopt.statespace import nns, greedy_tree
from ppopt.lp import lp_method
from ppopt.dp import dp_method
from ppopt.auxil import *

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
                partnum[i] = 0

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
                newcode = (i,partnum[i])  # determine which entry to record
                dic['parts'][newcode] = {'part_name': parts[i]}  # put the entry into dictionary
                well_parts.append(newcode)  # record part code in w
                reqvols[newcode] = part_vol  # record the required volume for the part - which is the same in BASIC
                partnum[i] += 1  # update number of parts of this class

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
            lp_method(w, fin, reord=None, caps=caps, maxtime=1)
        else:
            lp_method(w, fin, method[3:], caps=caps,maxtime=1)
    elif (method[0:2] == 'DP'):  # Dynamic Programming-based
        if (len(method) == 2):
            dp_method(w, fin, reord=None, caps=caps)
        else:
            dp_method(w, fin, method[3:], caps=caps)
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

    # PART 3.2 Define auxiliary variables
    added = np.zeros((len(w), len(w[0])))  # tells which parts were added to which well

    # PART 3.3 The first operation in fin
    added[fin[0].well][fin[0].part[0]] = 1
    fin[0].changed = True
    part_source = dic['parts'][fin[0].part]['part_name']
    part_vol = reqvols[fin[0].part]
    part_dest = [dic['constructs'][fin[0].well]['con_name']].copy()

    # PART 3.4 All later operations
    for i in range(1, len(fin)):
        # get operation cost
        cost = cost_func_with_w(fin[0:i], fin[i], w, added, caps)
        if (cost == 1):
            fin[i].changed = True
        added[fin[i].well][fin[i].part[0]] = 1

        # act accroding to cost
        if (cost == 0):  # if 0, tip is unchanged
            part_dest += [dic['constructs'][fin[i].well]['con_name']]
        else:  # if 1, there is new tip, so...
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


def moclo_execute(action_list, pipette):
    # Get tip to distribute DNA parts
    pipette.pick_up_tip()
    for action in action_list:
        # PART 1 Get directions
        source_well, destination_wells, volume, new_tip, air_gap = action
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
        # PART 5 Wash tip to recycle it
        p10_single.mix(2, 10, wash_0.bottom(0.5))
        p10_single.blow_out()
        p10_single.mix(2, 10, wash_1.bottom(0.5))
        p10_single.blow_out()

    # Drop the tip
    pipette.drop_tip()