# Estimation of automated protocol runtime by simulation
# WARNING: this code requires opentrons 4.7.0 to run
# This is NOT the version for which the rest of the package was written (3.21.1)
# To time the protocol, run 'opentrons_simulate -e estimate_runtimes.py'

from opentrons import protocol_api

# select which test to run
# - estimation of time to change a single tip by
# -- Distribute a DNA part across 96 wells with a fresh tip each time: WHICH_TEST='96 fresh'
# -- Distribute a DNA part across 96 wells, using one tip to serve two wells at a time: WHICH_TEST='96 reuse'
# - estimation of runtime for DNA part distribution in the test case from the DNA-BOT package
# -- Original DNA-BOT protocol, rewritten for Opentrons API 2.0: WHICH_TEST='DNA-BOT unoptimized'
# -- DNA-BOT protocol optimized by pipette_opt (DNA part order as read), rewritten for Opentrons API 2.0: WHICH_TEST='DNA-BOT optimized'
WHICH_TEST = 'DNA-BOT optimized'

# protocol metadata (for Opentrons API 2.0)
metadata = {
    'protocolName': 'My Protocol',
    'author': 'Name <email@address.com>',
    'description': 'Simple protocol to get started using OT2',
    'apiLevel': '2.11'
}

# dictionary of assembly master mixes to be made (original DNA-BOT code)
final_assembly_dict = {"A1": ["A7", "G7", "H7", "B8", "E8"], "B1": ["A7", "G7", "H7", "B8", "G8"],
                       "C1": ["A7", "G7", "H7", "A9", "E8"], "D1": ["A7", "G7", "H7", "A9", "G8"],
                       "E1": ["A7", "D9", "B8", "E9", "G9"], "F1": ["A7", "D9", "B8", "A10", "G9"],
                       "G1": ["A7", "D9", "A9", "E9", "G9"], "H1": ["A7", "D9", "A9", "A10", "G9"],
                       "A2": ["A7", "G7", "C10", "B8", "E8"], "B2": ["A7", "G7", "C10", "B8", "G8"],
                       "C2": ["A7", "G7", "C10", "A9", "E8"], "D2": ["A7", "G7", "C10", "A9", "G8"],
                       "E2": ["A7", "D9", "B8", "E9", "E10"], "F2": ["A7", "D9", "B8", "A10", "E10"],
                       "G2": ["A7", "D9", "A9", "E9", "E10"], "H2": ["B7", "D9", "A9", "A10", "E10"],
                       "A3": ["B7", "G10", "H10", "B8", "E8"], "B3": ["B7", "G10", "H10", "B8", "G8"],
                       "C3": ["B7", "G10", "H10", "A9", "E8"], "D3": ["B7", "G10", "H10", "A9", "G8"],
                       "E3": ["B7", "A11", "B8", "E9", "B11"], "F3": ["B7", "A11", "B8", "A10", "B11"],
                       "G3": ["B7", "A11", "A9", "E9", "B11"], "H3": ["B7", "A11", "A9", "A10", "B11"],
                       "A4": ["B7", "G10", "H7", "B8", "E8"], "B4": ["B7", "G10", "H7", "B8", "G8"],
                       "C4": ["B7", "G10", "H7", "A9", "E8"], "D4": ["B7", "G10", "H7", "A9", "G8"],
                       "E4": ["B7", "A11", "B8", "E9", "G9"], "F4": ["B7", "A11", "C8", "A10", "G9"],
                       "G4": ["C7", "A11", "A9", "E9", "G9"], "H4": ["C7", "A11", "B9", "A10", "G9"],
                       "A5": ["C7", "G10", "C10", "C8", "E8"], "B5": ["C7", "G10", "C10", "C8", "G8"],
                       "C5": ["C7", "G10", "C10", "B9", "E8"], "D5": ["C7", "G10", "C10", "B9", "G8"],
                       "E5": ["C7", "A11", "C8", "E9", "E10"], "F5": ["C7", "A11", "C8", "A10", "E10"],
                       "G5": ["C7", "A11", "B9", "E9", "E10"], "H5": ["C7", "A11", "B9", "A10", "E10"],
                       "A6": ["C7", "C11", "H10", "C8", "E8"], "B6": ["C7", "C11", "H10", "C8", "G8"],
                       "C6": ["C7", "C11", "H10", "B9", "E8"], "D6": ["C7", "C11", "H10", "B9", "G8"],
                       "E6": ["C7", "D11", "C8", "E9", "B11"], "F6": ["D7", "D11", "C8", "A10", "B11"],
                       "G6": ["D7", "D11", "B9", "E9", "B11"], "H6": ["D7", "D11", "B9", "A10", "B11"],
                       "A7": ["D7", "C11", "H7", "C8", "E8"], "B7": ["D7", "C11", "H7", "C8", "G8"],
                       "C7": ["D7", "C11", "H7", "B9", "E8"], "D7": ["D7", "C11", "H7", "B9", "G8"],
                       "E7": ["D7", "D11", "C8", "E9", "G9"], "F7": ["D7", "D11", "C8", "A10", "G9"],
                       "G7": ["D7", "D11", "B9", "E9", "G9"], "H7": ["D7", "D11", "B9", "A10", "G9"],
                       "A8": ["D7", "C11", "C10", "C8", "E8"], "B8": ["D7", "C11", "C10", "C8", "G8"],
                       "C8": ["D7", "C11", "C10", "B9", "F8"], "D8": ["D7", "C11", "C10", "B9", "H8"],
                       "E8": ["E7", "D11", "D8", "E9", "E10"], "F8": ["E7", "D11", "D8", "A10", "E10"],
                       "G8": ["E7", "D11", "C9", "F9", "E10"], "H8": ["E7", "D11", "C9", "B10", "E10"],
                       "A9": ["E7", "E11", "H10", "D8", "F8"], "B9": ["E7", "E11", "H10", "D8", "H8"],
                       "C9": ["E7", "E11", "H10", "C9", "F8"], "D9": ["E7", "E11", "H10", "C9", "H8"],
                       "E9": ["E7", "F11", "D8", "F9", "B11"], "F9": ["E7", "F11", "D8", "B10", "B11"],
                       "G9": ["E7", "F11", "C9", "F9", "B11"], "H9": ["E7", "F11", "C9", "B10", "B11"],
                       "A10": ["E7", "E11", "H7", "D8", "F8"], "B10": ["E7", "E11", "H7", "D8", "H8"],
                       "C10": ["E7", "E11", "H7", "C9", "F8"], "D10": ["F7", "E11", "A8", "C9", "H8"],
                       "E10": ["F7", "F11", "D8", "F9", "G9"], "F10": ["F7", "F11", "D8", "B10", "G9"],
                       "G10": ["F7", "F11", "C9", "F9", "G9"], "H10": ["F7", "F11", "C9", "B10", "H9"],
                       "A11": ["F7", "E11", "C10", "D8", "F8"], "B11": ["F7", "E11", "C10", "D8", "H8"],
                       "C11": ["F7", "E11", "C10", "C9", "F8"], "D11": ["F7", "E11", "D10", "C9", "H8"],
                       "E11": ["F7", "F11", "D8", "F9", "E10"], "F11": ["F7", "F11", "D8", "B10", "E10"],
                       "G11": ["F7", "F11", "C9", "F9", "E10"], "H11": ["F7", "F11", "C9", "B10", "F10"]}

# list of optimised pipette actions (DNA-BOT program optimised by pipette_opt)
action_list = (('A7', ['A1', 'B1', 'C1'], 1.5, 'once', 1), ('A7', ['E1', 'G1', 'F1', 'D1'], 1.5, 'once', 1),
               ('A7', ['A2', 'C2', 'B2', 'H1'], 1.5, 'once', 1), ('A7', ['E2', 'G2', 'F2', 'D2'], 1.5, 'once', 1),
               ('G7', ['A1', 'B1', 'D1', 'C1'], 1.5, 'once', 1), ('G7', ['B2', 'D2', 'C2', 'A2'], 1.5, 'once', 1),
               ('H7', ['A1', 'B1', 'C1'], 1.5, 'once', 1), ('H7', ['A4', 'C4', 'B4', 'D1'], 1.5, 'once', 1),
               ('H7', ['A7', 'C7', 'B7', 'D4'], 1.5, 'once', 1), ('H7', ['A10', 'C10', 'B10', 'D7'], 1.5, 'once', 1),
               ('B8', ['B2', 'A2', 'A1', 'B1'], 1.5, 'once', 1), ('B8', ['B3', 'A3', 'B4', 'A4'], 1.5, 'once', 1),
               ('E8', ['C2', 'C1', 'A1'], 1.5, 'once', 1), ('E8', ['C4', 'A7', 'C7', 'A4'], 1.5, 'once', 1),
               ('E8', ['A5', 'C3', 'A3', 'A2'], 1.5, 'once', 1), ('E8', ['A6', 'A8', 'C6', 'C5'], 1.5, 'once', 1),
               ('G8', ['D2', 'D1', 'B1'], 1.5, 'once', 1), ('G8', ['D4', 'B7', 'D7', 'B4'], 1.5, 'once', 1),
               ('G8', ['B5', 'D3', 'B3', 'B2'], 1.5, 'once', 1), ('G8', ['B6', 'B8', 'D6', 'D5'], 1.5, 'once', 1),
               ('A9', ['C2'], 1.5, 'once', 1), ('A9', ['D2'], 1.5, 'once', 1),
               ('A9', ['C3', 'C4', 'C1'], 1.5, 'once', 1), ('A9', ['D3', 'D4', 'D1'], 1.5, 'once', 1),
               ('D9', ['G2', 'F2', 'E2', 'H1'], 1.5, 'once', 1), ('D9', ['H2', 'E1', 'F1', 'G1'], 1.5, 'once', 1),
               ('B8', ['E1', 'F1', 'E2'], 1.5, 'once', 1), ('B8', ['E3', 'E4', 'F3', 'F2'], 1.5, 'once', 1),
               ('E9', ['G2', 'G1', 'E1'], 1.5, 'once', 1), ('E9', ['G3', 'E4', 'E3', 'E2'], 1.5, 'once', 1),
               ('E9', ['E5', 'E6', 'G5', 'G4'], 1.5, 'once', 1), ('E9', ['E7', 'E8', 'G7', 'G6'], 1.5, 'once', 1),
               ('G9', ['E7', 'G7', 'G4', 'G1'], 1.5, 'once', 1), ('G9', ['H7', 'F4', 'H1', 'E1'], 1.5, 'once', 1),
               ('G9', ['F10', 'F7', 'F1'], 1.5, 'once', 1), ('G9', ['G10', 'E10', 'H4', 'E4'], 1.5, 'once', 1),
               ('A10', ['H2', 'F2', 'F1'], 1.5, 'once', 1), ('A10', ['F6', 'H5', 'F5', 'F3'], 1.5, 'once', 1),
               ('A10', ['F7', 'F4', 'H4', 'H7'], 1.5, 'once', 1), ('A10', ['F8', 'H6', 'H3', 'H1'], 1.5, 'once', 1),
               ('A9', ['G2'], 1.5, 'once', 1), ('A9', ['G3', 'G4', 'G1'], 1.5, 'once', 1),
               ('A9', ['H3', 'H2', 'H1'], 1.5, 'once', 1), ('C10', ['B5', 'B8', 'D5', 'B2'], 1.5, 'once', 1),
               ('C10', ['C5', 'A8', 'A5', 'A2'], 1.5, 'once', 1), ('C10', ['C8', 'D8', 'A11', 'C2'], 1.5, 'once', 1),
               ('C10', ['B11', 'C11', 'D2'], 1.5, 'once', 1), ('E10', ['H5', 'F8', 'H2'], 1.5, 'once', 1),
               ('E10', ['G8', 'E5', 'G5', 'G2'], 1.5, 'once', 1), ('E10', ['H8', 'E11', 'F5', 'F2'], 1.5, 'once', 1),
               ('E10', ['F11', 'G11', 'E8', 'E2'], 1.5, 'once', 1), ('B7', ['A3', 'A4'], 1.5, 'once', 1),
               ('B7', ['B3', 'B4'], 1.5, 'once', 1), ('B7', ['C3', 'C4'], 1.5, 'once', 1),
               ('B7', ['D3', 'D4'], 1.5, 'once', 1), ('B7', ['E3', 'E4'], 1.5, 'once', 1),
               ('B7', ['F3'], 1.5, 'once', 1), ('B7', ['G3'], 1.5, 'once', 1), ('B7', ['H3', 'H2'], 1.5, 'once', 1),
               ('B7', ['F4'], 1.5, 'once', 1), ('G10', ['A3', 'A4'], 1.5, 'once', 1),
               ('G10', ['B3', 'B4'], 1.5, 'once', 1), ('G10', ['C3', 'C4'], 1.5, 'once', 1),
               ('G10', ['D3', 'D4'], 1.5, 'once', 1), ('G10', ['A5', 'C5'], 1.5, 'once', 1),
               ('G10', ['D5', 'B5'], 1.5, 'once', 1), ('H10', ['B6', 'D6', 'D3'], 1.5, 'once', 1),
               ('H10', ['A9', 'A6', 'C3'], 1.5, 'once', 1), ('H10', ['B9', 'B3'], 1.5, 'once', 1),
               ('H10', ['C9', 'D9', 'C6', 'A3'], 1.5, 'once', 1), ('A11', ['E3', 'E4'], 1.5, 'once', 1),
               ('A11', ['F3'], 1.5, 'once', 1), ('A11', ['G3'], 1.5, 'once', 1), ('A11', ['H3'], 1.5, 'once', 1),
               ('A11', ['G4'], 1.5, 'once', 1), ('A11', ['H4', 'F4'], 1.5, 'once', 1),
               ('A11', ['E5', 'G5'], 1.5, 'once', 1), ('A11', ['H5', 'F5'], 1.5, 'once', 1),
               ('B11', ['F6', 'H6', 'H3'], 1.5, 'once', 1), ('B11', ['E9', 'E6', 'G3'], 1.5, 'once', 1),
               ('B11', ['F9', 'F3'], 1.5, 'once', 1), ('B11', ['G9', 'H9', 'G6', 'E3'], 1.5, 'once', 1),
               ('C8', ['E5'], 1.5, 'once', 1), ('C8', ['F5'], 1.5, 'once', 1), ('C8', ['E6'], 1.5, 'once', 1),
               ('C8', ['F6'], 1.5, 'once', 1), ('C8', ['E7'], 1.5, 'once', 1), ('C8', ['F7', 'F4'], 1.5, 'once', 1),
               ('C7', ['G4'], 1.5, 'once', 1), ('C7', ['H4'], 1.5, 'once', 1), ('C7', ['A5', 'C5'], 1.5, 'once', 1),
               ('C7', ['B5', 'D5'], 1.5, 'once', 1), ('C7', ['G5', 'E5'], 1.5, 'once', 1),
               ('C7', ['H5', 'F5'], 1.5, 'once', 1), ('C7', ['B6', 'D6'], 1.5, 'once', 1),
               ('C7', ['C6', 'A6'], 1.5, 'once', 1), ('C7', ['E6'], 1.5, 'once', 1), ('B9', ['G5'], 1.5, 'once', 1),
               ('B9', ['H5'], 1.5, 'once', 1), ('B9', ['G6'], 1.5, 'once', 1), ('B9', ['H6'], 1.5, 'once', 1),
               ('B9', ['G7'], 1.5, 'once', 1), ('B9', ['H7', 'H4'], 1.5, 'once', 1), ('C8', ['A6'], 1.5, 'once', 1),
               ('C8', ['B6'], 1.5, 'once', 1), ('C8', ['A7'], 1.5, 'once', 1), ('C8', ['B7'], 1.5, 'once', 1),
               ('C8', ['A8', 'A5'], 1.5, 'once', 1), ('C8', ['B8', 'B5'], 1.5, 'once', 1),
               ('B9', ['C6'], 1.5, 'once', 1), ('B9', ['D6'], 1.5, 'once', 1), ('B9', ['C7'], 1.5, 'once', 1),
               ('B9', ['D7'], 1.5, 'once', 1), ('B9', ['C8', 'C5'], 1.5, 'once', 1),
               ('B9', ['D8', 'D5'], 1.5, 'once', 1), ('C11', ['A6'], 1.5, 'once', 1), ('C11', ['B6'], 1.5, 'once', 1),
               ('C11', ['C6'], 1.5, 'once', 1), ('C11', ['D6'], 1.5, 'once', 1), ('C11', ['A7'], 1.5, 'once', 1),
               ('C11', ['B7'], 1.5, 'once', 1), ('C11', ['C7'], 1.5, 'once', 1), ('C11', ['D7'], 1.5, 'once', 1),
               ('C11', ['A8'], 1.5, 'once', 1), ('C11', ['B8'], 1.5, 'once', 1), ('C11', ['C8', 'D8'], 1.5, 'once', 1),
               ('D11', ['E6'], 1.5, 'once', 1), ('D11', ['F6'], 1.5, 'once', 1), ('D11', ['G6'], 1.5, 'once', 1),
               ('D11', ['H6'], 1.5, 'once', 1), ('D11', ['E7'], 1.5, 'once', 1), ('D11', ['F7'], 1.5, 'once', 1),
               ('D11', ['G7'], 1.5, 'once', 1), ('D11', ['H7'], 1.5, 'once', 1), ('D11', ['G8', 'E8'], 1.5, 'once', 1),
               ('D11', ['H8', 'F8'], 1.5, 'once', 1), ('D7', ['F6'], 1.5, 'once', 1), ('D7', ['G6'], 1.5, 'once', 1),
               ('D7', ['H6'], 1.5, 'once', 1), ('D7', ['A7'], 1.5, 'once', 1), ('D7', ['B7'], 1.5, 'once', 1),
               ('D7', ['C7'], 1.5, 'once', 1), ('D7', ['D7'], 1.5, 'once', 1), ('D7', ['E7'], 1.5, 'once', 1),
               ('D7', ['F7'], 1.5, 'once', 1), ('D7', ['G7'], 1.5, 'once', 1), ('D7', ['H7'], 1.5, 'once', 1),
               ('D7', ['A8'], 1.5, 'once', 1), ('D7', ['B8'], 1.5, 'once', 1), ('D7', ['C8', 'D8'], 1.5, 'once', 1),
               ('F8', ['C9', 'A9'], 1.5, 'once', 1), ('F8', ['C10', 'A10'], 1.5, 'once', 1),
               ('F8', ['C11', 'A11', 'C8'], 1.5, 'once', 1), ('H8', ['D9', 'B9'], 1.5, 'once', 1),
               ('H8', ['D10', 'B11', 'D8'], 1.5, 'once', 1), ('H8', ['D11', 'B10'], 1.5, 'once', 1),
               ('E7', ['E8'], 1.5, 'once', 1), ('E7', ['H8', 'G8', 'F8'], 1.5, 'once', 1),
               ('E7', ['A9', 'C9'], 1.5, 'once', 1), ('E7', ['D9', 'B9'], 1.5, 'once', 1),
               ('E7', ['E9', 'F9', 'G9', 'H9'], 1.5, 'once', 1), ('E7', ['B10'], 1.5, 'once', 1),
               ('E7', ['C10', 'A10'], 1.5, 'once', 1), ('D8', ['E8'], 1.5, 'once', 1),
               ('D8', ['E9', 'F9'], 1.5, 'once', 1), ('D8', ['F10', 'E10'], 1.5, 'once', 1),
               ('D8', ['F11', 'E11', 'F8'], 1.5, 'once', 1), ('C9', ['H9', 'G9'], 1.5, 'once', 1),
               ('C9', ['H10', 'G10'], 1.5, 'once', 1), ('C9', ['H11', 'G11', 'H8', 'G8'], 1.5, 'once', 1),
               ('F9', ['E9'], 1.5, 'once', 1), ('F9', ['G9'], 1.5, 'once', 1), ('F9', ['E10'], 1.5, 'once', 1),
               ('F9', ['G10'], 1.5, 'once', 1), ('F9', ['E11'], 1.5, 'once', 1), ('F9', ['G11', 'G8'], 1.5, 'once', 1),
               ('B10', ['F9'], 1.5, 'once', 1), ('B10', ['F10'], 1.5, 'once', 1),
               ('B10', ['H10', 'H8'], 1.5, 'once', 1), ('B10', ['F11'], 1.5, 'once', 1),
               ('B10', ['H11', 'H9'], 1.5, 'once', 1), ('E11', ['A9', 'C9'], 1.5, 'once', 1),
               ('E11', ['D9', 'B9'], 1.5, 'once', 1), ('E11', ['A10', 'C10'], 1.5, 'once', 1),
               ('E11', ['B10'], 1.5, 'once', 1), ('E11', ['D10', 'D11', 'B11'], 1.5, 'once', 1),
               ('E11', ['C11', 'A11'], 1.5, 'once', 1), ('D8', ['A9'], 1.5, 'once', 1), ('D8', ['B9'], 1.5, 'once', 1),
               ('D8', ['A10'], 1.5, 'once', 1), ('D8', ['B10'], 1.5, 'once', 1), ('D8', ['A11'], 1.5, 'once', 1),
               ('D8', ['B11'], 1.5, 'once', 1), ('C9', ['C9'], 1.5, 'once', 1), ('C9', ['C10'], 1.5, 'once', 1),
               ('C9', ['C11'], 1.5, 'once', 1), ('C9', ['D11', 'D10', 'D9'], 1.5, 'once', 1),
               ('F11', ['E9'], 1.5, 'once', 1), ('F11', ['F9'], 1.5, 'once', 1), ('F11', ['G9'], 1.5, 'once', 1),
               ('F11', ['E10'], 1.5, 'once', 1), ('F11', ['F10'], 1.5, 'once', 1), ('F11', ['G10'], 1.5, 'once', 1),
               ('F11', ['E11'], 1.5, 'once', 1), ('F11', ['F11'], 1.5, 'once', 1), ('F11', ['G11'], 1.5, 'once', 1),
               ('F11', ['H11', 'H10', 'H9'], 1.5, 'once', 1), ('F7', ['D10', 'D11'], 1.5, 'once', 1),
               ('F7', ['E10'], 1.5, 'once', 1), ('F7', ['F10'], 1.5, 'once', 1), ('F7', ['G10'], 1.5, 'once', 1),
               ('F7', ['A11'], 1.5, 'once', 1), ('F7', ['B11'], 1.5, 'once', 1), ('F7', ['C11'], 1.5, 'once', 1),
               ('F7', ['E11'], 1.5, 'once', 1), ('F7', ['F11'], 1.5, 'once', 1), ('F7', ['G11'], 1.5, 'once', 1),
               ('F7', ['H11', 'H10'], 1.5, 'once', 1), ('A8', ['D10'], 1.5, 'once', 1), ('H9', ['H10'], 1.5, 'once', 1),
               ('D10', ['D11'], 1.5, 'once', 1), ('F10', ['H11'], 1.5, 'once', 1))


# Function for distributing DNA parts in optimised DNA-BOT program
def basic_execute(action_list, pipette, magbead_plate, destination_plate):
    for action in action_list:
        # PART 1 Get directions
        source_well, destination_wells, volume, new_tip, air_gap = action

        # PART 2 Get new tip
        pipette.pick_up_tip()

        # PART 3 Aspirate
        # 3a) with air gaps
        if (air_gap != 0):
            for i in range(0, len(destination_wells) - 1):
                pipette.aspirate(volume, magbead_plate.well(source_well))
                pipette.air_gap(air_gap)
            pipette.aspirate(volume, magbead_plate.well(source_well))
        # 3b) no air gaps
        else:
            pipette.aspirate(volume * len(destination_wells), magbead_plate.well(source_well))

        # PART 4 Distribute
        for i in range(0, len(destination_wells) - 1):
            pipette.dispense(volume + air_gap, destination_plate.well(destination_wells[i]))
        pipette.dispense(volume, destination_plate.well(destination_wells[- 1]))

        # PART 5 Drop used tip
        pipette.drop_tip()

    return


# ---OPENTRONS PROTOCOL---
def run(protocol: protocol_api.ProtocolContext):
    if (WHICH_TEST == '96 fresh' or WHICH_TEST == '96 reused'):
        # Define tip rack location
        tiprack = protocol.load_labware('opentrons_96_tiprack_10ul', '1')

        # Define DNA part source plate
        part_plate = protocol.load_labware('corning_96_wellplate_360ul_flat', '6')

        # Define well plate with constructs being prepared
        construct_plate = protocol.load_labware('corning_96_wellplate_360ul_flat', '7')

        # Define pipette
        left_pipette = protocol.load_instrument(
            'p10_single', 'left', tip_racks=[tiprack])

        if (WHICH_TEST == '96 fresh'):
            for i in range(0, 96):
                left_pipette.pick_up_tip()
                left_pipette.aspirate(1, part_plate['A1'])
                left_pipette.dispense(1, construct_plate.wells()[i])
                left_pipette.drop_tip()
        else:
            for i in range(0, 48):
                left_pipette.pick_up_tip()
                left_pipette.aspirate(2, part_plate['A1'])
                left_pipette.dispense(1, construct_plate.wells()[2 * i])
                left_pipette.dispense(1, construct_plate.wells()[2 * i + 1])
                left_pipette.drop_tip()
    else:
        # Define tip racks (5 racks needed to satisfy demand for up to 440 tips)
        tiprack1 = protocol.load_labware('opentrons_96_tiprack_10ul', '1')
        tiprack4 = protocol.load_labware('opentrons_96_tiprack_10ul', '2')
        tiprack5 = protocol.load_labware('opentrons_96_tiprack_10ul', '3')
        tiprack6 = protocol.load_labware('opentrons_96_tiprack_10ul', '4')
        tiprack7 = protocol.load_labware('opentrons_96_tiprack_10ul', '5')

        # Define well plate with magnetic beads (DNA part source) - made to be Corning, not custom-defined FrameStar, for simplicity
        magbead_plate = protocol.load_labware('corning_96_wellplate_360ul_flat', '6')

        # Define well plate with constructs being prepared- made to be Corning, not custom-defined FrameStar, for simplicity
        destination_plate = protocol.load_labware('corning_96_wellplate_360ul_flat', '7')

        # Define pipette
        left_pipette = protocol.load_instrument(
            'p10_single', 'left', tip_racks=[tiprack1, tiprack4, tiprack5, tiprack6, tiprack7])

        if (WHICH_TEST == 'DNA-BOT unoptimized'):
            for key, values in list(final_assembly_dict.items()):
                left_pipette.transfer(1.5, magbead_plate.wells(values[0], values[1], values[2], values[3], values[4]),
                                      destination_plate.wells(key), mix_after=(1, 3), new_tip='always')
        else:
            basic_execute(action_list, left_pipette, magbead_plate, destination_plate)

    return