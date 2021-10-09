# startstop - API for the DNA-BOT automated assmembly package (Start-Stop standard)
# v0.1.0, 31.5.21

"""
Make the following changes to the DNA-BOT files:
(Here, the indents are used to aid legibility; appropriate formatting after inserting is required)

In dnabot_app.py:
- After all import statements, insert:
    # Importing API for ppopt optimisation of pipette tip consumption
    from ppopt.basic import basic_part_transfer_actions

- Line 103 (after the 'calculate OT2 script variables' section). Insert (method selection the same as for Start-Stop):
    print('pipette_opt: optimising part distribution order...')
    ppopt_action_list = basic_part_transfer_actions(method=..., final_assembly_dict=final_assembly_dict, part_vol=1.5, pipette_vol = 10)

- Line 114 - add a kwarg action_list=ppopt_action_list to get:
    generate_ot2_script(F_ASSEMBLY_FNAME, os.path.join(
        template_dir_path, F_ASSEMBLY_TEMP_FNAME),
        final_assembly_dict=final_assembly_dict,
        tiprack_num=final_assembly_tipracks,action_list=ppopt_action_list)



In assembly_template.py:
- In final_assembly() function:
    Comment the whole 'part transfers' section
    After it, put  a line:
        basic_execute(action_list, pipette, magbead_plate, destination_plate)
"""

from ppopt.statespace import nns, greedy_tree
from ppopt.lp import lp_method
from ppopt.dp import dp_method
from ppopt.auxil import *
from ppopt.visualise import rec

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
                partnum[i] = 0

            # determine part name
            partloc = parts[i]

            # check if such location is already in the dictionary
            match = False
            for partcode in dic['parts'].keys():
                # note that the same part in two different positions in the consrtucts is considered two different parts
                if (partloc == dic['parts'][partcode]['part_liqloc'] and i == partcode[0]):  # yes, record the part code in w
                    well_parts.append(partcode)
                    match = True
                    break

            if (not match):  # if no, update the dictionary
                newcode = (i,partnum[i])  # determine which entry to record
                dic['parts'][newcode] = {'part_liqloc': parts[i]}  # put the entry into dictionary
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
            lp_method(w, fin, method[3:], caps=caps, maxtime=1)
    elif (method[0:2] == 'DP'):  # LP-based
        if (len(method) == 2):
            dp_method(w, fin, reord=None, caps=caps)
        else:
            dp_method(w, fin, reord=method[3:], caps=caps)
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