'''
Course: Introduction to Bioinformatics
Assignment: Midterm Milestone #2
Date: Nov 22, 2024

This python3 program is designed to recursively globally align a set
of protein sequences into one alignment.
'''

#Make sure protein sequences have valid characters
def ValidateProteins(proteinSeqs, nameOfSeqs, substitutionMatrix):
    validSeqs = substitutionMatrix.keys()
    
    # make sure only seqs in proteinSeqs are the seqs that exist
    for i, seq in enumerate(proteinSeqs):
        if any(sequences not in validSeqs for sequences in seq):
            if nameOfSeqs:
                print(f"The sequence " + str(nameOfSeqs[i]) + " is not valid. Please try again.")
                return False
            else:
                print(f"The sequence " + str(proteinSeqs[i]) + " is not valid. Please try again.")
                return False
    return True


#This is our team's substitution matrix we are using to determine alignment for protein sequences
# '-': {'-': 3.4253431744398948, 'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 0, 'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 0, 'S': 0, 'T': 0, 'V': 0, 'W': 0, 'X': 0, 'Y': 0},
substitutionMatrix = {
    '-': {'-': 3.4253431744398948, 'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 0, 'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 0, 'S': 0, 'T': 0, 'V': 0, 'W': 0, 'X': 0, 'Y': 0},
    'A': {'-': -2.119818318223286, 'A': 2.1713614795485054, 'C': -0.0130122468148803, 'D': -0.35985306763820957, 'E': -0.1762968101545147, 'F': -1.9633669357201928, 'G': 0.22747035077434463, 'H': -1.8446496107785948, 'I': -1.0001288311236913, 'K': -0.447400767811872, 'L': -1.7582812430898471, 'M': -1.288171341526304, 'N': -0.7618694495981007, 'P': -0.05737699731902333, 'Q': -0.034224209254968795, 'R': -0.5073105470793311, 'S': 0.34311276007850533, 'T': -0.27095814510528177, 'V': -0.30796733770868173, 'W': -2.35552714223542, 'X': -0.6696357326634829, 'Y': -1.9536253448909018},
    'C': {'-': -1.6419667605195312, 'A': 0, 'C': 4.787289030363239, 'D': -2.423983405057925, 'E': -1.1561147771080602, 'F': -1.1582298389124077, 'G': -0.5507278911690605, 'H': -2.8860200451984555, 'I': -0.11558850273865542, 'K': -0.9326701384972766, 'L': -0.08247771188778082, 'M': -0.8778342724171765, 'N': -0.09816742789878038, 'P': -0.9380134595715865, 'Q': -0.4991403833027363, 'R': -1.0158756111076142, 'S': 0.2048333852526584, 'T': 0.3105529304269452, 'V': 0.10404163985924456, 'W': -2.097729384767773, 'X': 2.5267614801400207, 'Y': 0.7171435065738474},
    'D': {'-': -1.257921399893065, 'A': 0, 'C': 0, 'D': 2.8802128249487686, 'E': 1.1021555271492218, 'F': -4.043879695888468, 'G': -0.48210745150007517, 'H': -1.0478472565954735, 'I': -3.0538424437534046, 'K': -0.6472332974691781, 'L': -2.8697238648028467, 'M': -2.395163129021424, 'N': 0.689673376737642, 'P': -0.16065272752355453, 'Q': 0.17629989494698706, 'R': -0.37512995258690573, 'S': -0.5247391167710956, 'T': -0.5400804699807965, 'V': -2.792165961329269, 'W': -2.3115086757995127, 'X': 2.312982189108281, 'Y': -1.1868071671451572},
    'E': {'-': -2.169263631689073, 'A': 0, 'C': 0, 'D': 0, 'E': 2.6494703585318207, 'F': -2.9893732486980933, 'G': -0.8240534399643951, 'H': -1.0233721128524147, 'I': -2.40286878056756, 'K': -0.4465752422987094, 'L': -2.3767886812182084, 'M': -2.2571034147436473, 'N': 0.34880644220611545, 'P': -0.1512390342587769, 'Q': 0.6020955145992908, 'R': -0.010091403853008774, 'S': -0.03117870609357994, 'T': -0.4675114448944255, 'V': -1.6416844178437868, 'W': -1.2010566374810614, 'X': -0.46192320663214537, 'Y': -1.4306076661659368},
    'F': {'-': -2.3689731350508736, 'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 3.627650665607347, 'G': -4.210152546186241, 'H': -2.091940032120465, 'I': -0.2200632261011195, 'K': -2.849222194159748, 'L': 0.25622554439310036, 'M': -0.6079134809981425, 'N': -2.6683914936550197, 'P': -2.940410314820683, 'Q': -2.226235733444881, 'R': -2.4265525994602455, 'S': -2.3847520274370186, 'T': -2.5553030571844872, 'V': -0.8340029194228896, 'W': -0.886225279574061, 'X': -0.01067265049854963, 'Y': 1.0360493907446329},
    'G': {'-': -1.0590970056477471, 'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 3.2437588286935366, 'H': -1.8582662564801993, 'I': -3.2880193585529804, 'K': -1.2510186798016134, 'L': -3.2656939370603006, 'M': -2.8428489144678446, 'N': 0.12713092849294444, 'P': -1.0698281969953853, 'Q': -0.8218050424747285, 'R': -0.7363759352952864, 'S': -0.11830693210299645, 'T': -0.8532906611894912, 'V': -1.4175795936458977, 'W': -3.1116536272993875, 'X': -0.17521875607685383, 'Y': -2.6099883204798817},
    'H': {'-': -2.1367038817748, 'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 3.847283368079261, 'I': -1.8162049242960856, 'K': -1.1099562816002404, 'L': -2.1506372968962006, 'M': -2.6982737666688985, 'N': 0.44377016746187775, 'P': -2.2930794870537747, 'Q': 0.3641215577544446, 'R': -0.5563309005293989, 'S': -1.3514397336753499, 'T': -1.5294115240524637, 'V': -1.4583284259909812, 'W': -2.0598495010966835, 'X': 3.546939362077651, 'Y': 0.3496205938030878},
    'I': {'-': -2.5913427557213176, 'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 'I': 3.086782049262499, 'K': -1.4765179499197123, 'L': 0.8729875091771264, 'M': 1.0167003132909547, 'N': -1.9729526539713222, 'P': -2.300782727736586, 'Q': -1.425733985378439, 'R': -1.134719296386524, 'S': -1.4106732749982431, 'T': -1.25664007420692, 'V': 1.4157028994338177, 'W': -2.0937163222611495, 'X': 0.2088464477592817, 'Y': -0.2820176787146727},
    'K': {'-': -2.107745485922711, 'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 'I': 0, 'K': 2.563490592765242, 'L': -1.6404691572496923, 'M': -1.5596583123909213, 'N': -0.1306913548953145, 'P': -0.6736526029012139, 'Q': 0.624808307022226, 'R': 1.0531410180194074, 'S': -0.6439685002141643, 'T': -0.6482970044080562, 'V': -1.8251236579690227, 'W': -2.3858895762552415, 'X': 0.9273996003582483, 'Y': -1.3975466860920285},
    'L': {'-': -2.024245757186028, 'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 2.4097442570223904, 'M': 0.7151743655022903, 'N': -1.4476178735964418, 'P': -2.6115095254200487, 'Q': -1.2524213830729047, 'R': -1.511409356448346, 'S': -1.6591632662002216, 'T': -1.530944973240256, 'V': 0.07011035525168247, 'W': -1.4615282309765183, 'X': -1.01584051885749, 'Y': -0.49720277756429},
    'M': {'-': -1.6812344610667995, 'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 0, 'M': 3.965048120753108, 'N': -1.5241182205113628, 'P': -1.358550755997805, 'Q': -1.8989191195621575, 'R': -0.45717146770062966, 'S': -1.1236619837385786, 'T': -0.4899365236122178, 'V': 0.1928967555941452, 'W': -0.18781451400009044, 'X': 1.006688510162888, 'Y': 0.5780191452962775},
    'N': {'-': -1.4738863618645017, 'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 0, 'M': 0, 'N': 2.8627765548187636, 'P': -2.007727210283283, 'Q': 0.09459301949916628, 'R': -0.08154252875660851, 'S': 0.4297226660795161, 'T': 0.5216321094850765, 'V': -1.6241165301503464, 'W': -2.6748297677159583, 'X': 0.5960241425771345, 'Y': -0.42676772229610777},
    'P': {'-': -0.24910927185411424, 'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 0, 'M': 0, 'N': 0, 'P': 3.8736654791309997, 'Q': -0.6530384778111278, 'R': -0.8748093204221488, 'S': -0.11976389641068835, 'T': -0.9997175697372999, 'V': -1.3180214512397195, 'W': -2.4223266664311196, 'X': 0.6172016977555178, 'Y': -3.9997711978682595},
    'Q': {'-': -1.5904082284020455, 'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 0, 'M': 0, 'N': 0, 'P': 0, 'Q': 2.9698983834239234, 'R': 0.5773407346440468, 'S': -0.395675779072049, 'T': -0.5988608282474596, 'V': -1.2444188670127694, 'W': -3.7727240863514018, 'X': 2.147222662082563, 'Y': -0.8603346890955187},
    'R': {'-': -1.2319535116266533, 'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 0, 'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 3.0553016864183427, 'S': -0.16944751463512184, 'T': 0.07641781611259255, 'V': -1.6434879680399346, 'W': -0.9882897794357143, 'X': 0.5393395462194903, 'Y': -0.8592098302707846},
    'S': {'-': -0.9640222738252485, 'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 0, 'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 0, 'S': 2.4792888615366744, 'T': 0.7391949574772038, 'V': -1.0438063587744237, 'W': -2.8663406055504446, 'X': 1.5588414511339426, 'Y': -1.4487159000441394},
    'T': {'-': -1.4480458203026731, 'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 0, 'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 0, 'S': 0, 'T': 2.959874210657771, 'V': -0.09993346387042548, 'W': -2.7434459165680254, 'X': 1.6817361401163615, 'Y': -1.143823377318833},
    'V': {'-': -2.235531636624174, 'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 0, 'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 0, 'S': 0, 'T': 0, 'V': 2.531102710524824, 'W': -2.484648503620702, 'X': -0.36265797924209175, 'Y': -1.4226827503143094},
    'W': {'-': -2.9749948112670133, 'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 0, 'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 0, 'S': 0, 'T': 0, 'V': 0, 'W': 5.689996730540899, 'X': 1.902270615232227, 'Y': 1.449044146121638},
    'X': {'-': -0.538130949534991, 'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 0, 'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 0, 'S': 0, 'T': 0, 'V': 0, 'W': 0, 'X': 10.696686481582333, 'Y': 1.3792738678174639},
    'Y': {'-': -1.4892213490540445, 'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 0, 'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 0, 'S': 0, 'T': 0, 'V': 0, 'W': 0, 'X': 0, 'Y': 4.442322319141567}
}

gapPenalty = -2

# BLOSUM62 = {
#     '-': {'-': 3.4253431744398948, 'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 0, 'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 0, 'S': 0, 'T': 0, 'V': 0, 'W': 0, 'X': 0, 'Y': 0},
#     'A': {'A': 4, 'C': 0, 'D': -2, 'E': -1, 'F': -2, 'G': 0, 'H': -1, 'I': -1, 'K': -1, 'L': -1, 'M': -1, 'N': -2, 'P': -1, 'Q': -1, 'R': -1, 'S': 1, 'T': 0, 'V': 0, 'W': -3, 'Y': -2},
#     'C': {'A': 0, 'C': 9, 'D': -3, 'E': -4, 'F': -2, 'G': -3, 'H': -3, 'I': -1, 'K': -3, 'L': -1, 'M': -1, 'N': -3, 'P': -2, 'Q': -3, 'R': -3, 'S': -1, 'T': -1, 'V': -1, 'W': -2, 'Y': -2},
#     'D': {'A': -2, 'C': -3, 'D': 6, 'E': 2, 'F': -3, 'G': -1, 'H': 0, 'I': -3, 'K': -1, 'L': -3, 'M': -2, 'N': 1, 'P': -1, 'Q': 0, 'R': 0, 'S': 0, 'T': -1, 'V': -2, 'W': -4, 'Y': -3},
#     'E': {'A': -1, 'C': -4, 'D': 2, 'E': 5, 'F': -3, 'G': -2, 'H': 0, 'I': -3, 'K': 1, 'L': -2, 'M': -2, 'N': 0, 'P': -1, 'Q': 2, 'R': 1, 'S': 0, 'T': -1, 'V': -2, 'W': -3, 'Y': -2},
#     'F': {'A': -2, 'C': -2, 'D': -3, 'E': -3, 'F': 6, 'G': -3, 'H': -1, 'I': 0, 'K': -2, 'L': 0, 'M': 0, 'N': -3, 'P': -2, 'Q': -2, 'R': -2, 'S': -1, 'T': -1, 'V': 0, 'W': 1, 'Y': 3},
#     'G': {'A': 0, 'C': -3, 'D': -1, 'E': -2, 'F': -3, 'G': 6, 'H': -2, 'I': -4, 'K': -2, 'L': -4, 'M': -3, 'N': -1, 'P': -2, 'Q': -2, 'R': -2, 'S': 0, 'T': -2, 'V': -3, 'W': -2, 'Y': -3},
#     'H': {'A': -1, 'C': -3, 'D': 0, 'E': 0, 'F': -1, 'G': -2, 'H': 8, 'I': -2, 'K': -1, 'L': -2, 'M': -2, 'N': 1, 'P': -1, 'Q': 0, 'R': 0, 'S': -1, 'T': -2, 'V': -2, 'W': -1, 'Y': 2},
#     'I': {'A': -1, 'C': -1, 'D': -3, 'E': -3, 'F': 0, 'G': -4, 'H': -2, 'I': 4, 'K': -1, 'L': 2, 'M': 1, 'N': -3, 'P': -2, 'Q': -1, 'R': -1, 'S': -2, 'T': -1, 'V': 3, 'W': -3, 'Y': -1},
#     'K': {'A': -1, 'C': -3, 'D': -1, 'E': 1, 'F': -2, 'G': -2, 'H': -1, 'I': -1, 'K': 5, 'L': -1, 'M': -1, 'N': 0, 'P': -1, 'Q': 1, 'R': 2, 'S': 0, 'T': -1, 'V': -2, 'W': -3, 'Y': -2},
#     'L': {'A': -1, 'C': -1, 'D': -3, 'E': -2, 'F': 0, 'G': -4, 'H': -2, 'I': 2, 'K': -1, 'L': 4, 'M': 2, 'N': -3, 'P': -2, 'Q': -1, 'R': -1, 'S': -2, 'T': -1, 'V': 1, 'W': -2, 'Y': -1},
#     'M': {'A': -1, 'C': -1, 'D': -2, 'E': -2, 'F': 0, 'G': -3, 'H': -2, 'I': 1, 'K': -1, 'L': 2, 'M': 5, 'N': -2, 'P': -2, 'Q': 0, 'R': 0, 'S': -1, 'T': -1, 'V': 1, 'W': -1, 'Y': -1},
#     'N': {'A': -2, 'C': -3, 'D': 1, 'E': 0, 'F': -3, 'G': -1, 'H': 1, 'I': -3, 'K': 0, 'L': -3, 'M': -2, 'N': 6, 'P': -1, 'Q': 0, 'R': 0, 'S': 0, 'T': -1, 'V': -3, 'W': -4, 'Y': -3},
#     'P': {'A': -1, 'C': -2, 'D': -1, 'E': -1, 'F': -2, 'G': -2, 'H': -1, 'I': -2, 'K': -1, 'L': -2, 'M': -2, 'N': -1, 'P': 7, 'Q': -1, 'R': -1, 'S': -1, 'T': -1, 'V': -2, 'W': -3, 'Y': -2},
#     'Q': {'A': -1, 'C': -3, 'D': 0, 'E': 2, 'F': -2, 'G': -2, 'H': 0, 'I': -1, 'K': 1, 'L': -1, 'M': 0, 'N': 0, 'P': -1, 'Q': 5, 'R': 1, 'S': 0, 'T': -1, 'V': -2, 'W': -3, 'Y': -2},
#     'R': {'A': -1, 'C': -3, 'D': 0, 'E': 1, 'F': -2, 'G': -2, 'H': 0, 'I': -1, 'K': 2, 'L': -1, 'M': 0, 'N': 0, 'P': -1, 'Q': 1, 'R': 5, 'S': 0, 'T': -1, 'V': -2, 'W': -3, 'Y': -2},
#     'S': {'A': 1, 'C': -1, 'D': 0, 'E': 0, 'F': -1, 'G': 0, 'H': -1, 'I': -2, 'K': 0, 'L': -2, 'M': -1, 'N': 0, 'P': -1, 'Q': 0, 'R': 0, 'S': 4, 'T': 1, 'V': -2, 'W': -3, 'Y': -2},
#     'T': {'A': 0, 'C': -1, 'D': -1, 'E': -1, 'F': -1, 'G': -2, 'H': -2, 'I': -1, 'K': -1, 'L': -1, 'M': -1, 'N': -1, 'P': -1, 'Q': -1, 'R': -1, 'S': 1, 'T': 5, 'V': 1, 'W': -2, 'Y': -2},
#     'V': {'A': 0, 'C': -1, 'D': -2, 'E': -2, 'F': 0, 'G': -3, 'H': -2, 'I': 3, 'K': -2, 'L': 1, 'M': 1, 'N': -3, 'P': -2, 'Q': -2, 'R': -2, 'S': -2, 'T': 1, 'V': 4, 'W': -3, 'Y': -1},
#     'W': {'A': -3, 'C': -2, 'D': -4, 'E': -3, 'F': 1, 'G': -2, 'H': -1, 'I': -3, 'K': -3, 'L': -2, 'M': -1, 'N': -4, 'P': -3, 'Q': -3, 'R': -3, 'S': -3, 'T': -2, 'V': -3, 'W': 11, 'Y': 2},
#     'Y': {'A': -2, 'C': -2, 'D': -3, 'E': -2, 'F': 3, 'G': -3, 'H': 2, 'I': -1, 'K': -2, 'L': -1, 'M': -1, 'N': -3, 'P': -2, 'Q': -2, 'R': -2, 'S': -2, 'T': -2, 'V': -1, 'W': 2, 'Y': 7}
# }

def GloballyAlign(protein1, protein2):
    #Initialize matrix based on protein1 length and protein2 length
    #Make it one longer to include a gap at beginning of sequences in matrix
    protein1 = "-" + protein1
    protein2 = "-" + protein2
    matrixWidth = len(protein1)
    matrixHeight = len(protein2)
    matrix = [[0 for x in range(matrixHeight + 1)] for y in range(matrixWidth + 1)]

    optimalPathSeq1 = [] # character in protein 1 (i)
    optimalPathSeq2 = [] # character in protein 2 (j)

    # Edit cells in first row and column based on gap penalty
    for i in range(0, matrixWidth):
        #matrix at [0][0] should be 0
        if i == 0:
            matrix[i][0] = 0
        else:
            # matrix[i][0] = matrix[i - 1][0] + substitutionMatrix[protein1[i - 1]].get('-')
            # there is no gap penalty in our matrix so I will introduce a hard-coded gap penalty
            matrix[i][0] = matrix[i - 1][0] + gapPenalty

    for j in range(0, matrixHeight):
        if j == 0:
            matrix[0][j] = 0
        else:
            matrix[0][j] = matrix[0][j - 1] + gapPenalty
    
    # Edit the rest of the matrix cells
    for i in range(1, matrixWidth):
        for j in range(1, matrixHeight):            
            # For each cell (i, j), calculate:
            # (a) Cell (i – 1, j) – Gap penalty
            a = matrix[i - 1][j] + gapPenalty
            
            # (b) Cell (i, j – 1) – Gap penalty
            b = matrix[i][j - 1] + gapPenalty

            # (c) Cell (i – 1, j – 1) ± Match/mismatch
            # find match or mismatch
            #if it's an x, add gap penalty
            #DIAGONAL
            if protein1[i] == 'x' or protein2[j] == 'x':
                c = matrix[i - 1][j - 1] + gapPenalty
            else:
                #match or mismatch
                c = matrix[i - 1][j - 1] + substitutionMatrix[protein1[i]].get(protein2[j])

            # Then fill in the largest of these values
            matrix[i][j] = max(a, b, c)

    # backward pass time!!
    # final cell is as long as the shortest protein length
    i = matrixWidth - 1
    j = matrixHeight - 1
    current = matrix[i][j]

    optimalPathSeq1 = [protein1[i]] # character in protein 1 (i)
    optimalPathSeq2 = [protein2[j]] # character in protein 2 (j)

    while i > 0 and j > 0:
        top = matrix[i][j-1]
        left = matrix[i-1][j]
        diag = matrix[i-1][j-1]
        
        # chose best optimal path
        isDiag = False
        isLeft = False
        isTop = False
        #diagonal
        #deal with x or gaps first
        if (protein1[i] == 'x' or protein2[j] == 'x'):
            if(protein1[i] == 'x'):
                #add gap on left
                isTop = True
            else:
                isLeft = True
        #check if diagonal leads to current cell
        elif (diag + substitutionMatrix[protein1[i]].get(protein2[j]) == current):
            # possiblePaths.append((protein1[i], protein2[j]))
            isDiag = True
        #now check left
        elif (left + gapPenalty):
            isLeft = True
        #check top
        elif (top + gapPenalty):
            isTop = True
        
        #if diagonal, take the diagonal path
        if isDiag:
            optimalPathSeq1.append(protein1[i-1])
            optimalPathSeq2.append(protein2[j-1])
            current = matrix[i-1][j-1]
        #else if left, take the left path
        elif isLeft:
            optimalPathSeq1.append(protein1[i-1])
            optimalPathSeq2.append(protein2[j])
            current = matrix[i-1][j]
        #else if top, take the top path
        elif isTop:
            optimalPathSeq1.append(protein1[i])
            optimalPathSeq2.append(protein2[j-1])
            current = matrix[i-1][j-1]
        
        i -= 1
        j -= 1

    # add gaps for rest of sequence
    while i > 0:
        optimalPathSeq1.append(protein1[i-1])
        optimalPathSeq2.append('-')
        i -= 1
    
    while j > 0:
        optimalPathSeq1.append('-')
        optimalPathSeq2.append(protein2[j-1])
        j -= 1

    #remove gaps at front
    while optimalPathSeq1[len(optimalPathSeq1) - 1] == '-' and optimalPathSeq2[len(optimalPathSeq2) - 1] == '-':
        optimalPathSeq1.pop()
        optimalPathSeq2.pop()
    
    optimalPathSeq1 = ''.join(reversed(optimalPathSeq1))
    optimalPathSeq2 = ''.join(reversed(optimalPathSeq2))
    
    return optimalPathSeq1, optimalPathSeq2


def AlignSequences(proteinSeqs):
    # Analyze a set of protein sequences to find which pair shows highest similarity
    # Return the two sequences with the most matching protein sequences
    similarity = 0
    tempSimilarity = 0
    highestMatchIndex1 = -1
    highestMatchIndex2 = -1

    for i in range(len(proteinSeqs)):
        for j in range(i + 1, len(proteinSeqs)):
            protein1 = proteinSeqs[i]
            protein2 = proteinSeqs[j]
            if len(protein1) < len(protein2):
                smallerProteinLen = len(protein1)
            else:
                smallerProteinLen = len(protein2)

            similarity = 0  # Reset similarity
            for charIndex in range(smallerProteinLen):
                if protein1[charIndex] == protein2[charIndex]:
                    similarity += 1

            if similarity > tempSimilarity:
                highestMatchIndex1 = i
                highestMatchIndex2 = j
                tempSimilarity = similarity  # Store highest similarity so far in tempSimilarity

    # globally align those two sequences using your team’s substitution matrix
    aligned1, aligned2 = GloballyAlign(protein1, protein2)
    finalAligned = AlignTwoSequencesIntoOne(aligned1, aligned2, proteinSeqs)

    # Replace the two protein sequences with their alignment in the original set of sequences
    newProteinSeqs = []
    alignmentPlaced = False
    for i in range(len(proteinSeqs)):
        if not alignmentPlaced and (i == highestMatchIndex1 or i == highestMatchIndex2):
            #Append singular final aligned protein sequence to new protein sequences list
            newProteinSeqs.append(finalAligned)
            alignmentPlaced = True
        elif i != highestMatchIndex1 and i != highestMatchIndex2:
            newProteinSeqs.append(proteinSeqs[i])
    proteinSeqs = newProteinSeqs

    return newProteinSeqs

def AlignTwoSequencesIntoOne(protein1, protein2, proteinSeqs):
    aligned = ""

    shortestProtein = protein2
    if protein1 >= protein2:
        shortestProtein = protein1

    for i, aa in enumerate(shortestProtein):
        if protein1[i] == protein2[i]:
            aligned += aa 
        else:
            aminoAcidsPresent = {}
            for sequence in proteinSeqs:
                if i < len(sequence):
                    aa = sequence[i]
                    if aa in aminoAcidsPresent:
                        aminoAcidsPresent[aa] += 1
                    else:
                        aminoAcidsPresent[aa] = 1
            if aminoAcidsPresent:
                highestAACount = max(aminoAcidsPresent, key=aminoAcidsPresent.get)
            else:
                highestAACount = "-"
            aligned += highestAACount

    return aligned


def main():
    #Read input of FASTA formatted protein sequences
    inputFile = open("alignment.txt", "r").readlines()

    #Store the protein sequences along with the names of those sequences
    nameOfSeqs = []
    proteinSeqs = []
    readNextLine = False
    isFasta = False

    #Sort through input file
    for eachLine in inputFile:
        if eachLine.startswith(">"):
            isFasta = True
        if isFasta:
            #Save name for the file
            if readNextLine:
                #Save dna seq associated with that line and find transversions/transitions
                proteinSeqs.append(eachLine.strip())
                readNextLine = False
            if ">" in eachLine:
                #Skip the first char and save name to array if the bounds allow it
                nameOfSeqs.append(eachLine[1 : len(eachLine)].strip())
                #read the next line and store the DNA sequence
                readNextLine = True
        else:
            if eachLine.strip():
                proteinSeqs.append(eachLine.strip())
                
    #Validate protein sequences
    validProteinSeq = ValidateProteins(proteinSeqs, nameOfSeqs, substitutionMatrix)
    if not validProteinSeq:
        return

    #re-iterate the steps until all sequences have been included in a single multiple alignment.
    loopFound = False
    while len(proteinSeqs) > 1 and not loopFound:
        proteinSeqs = AlignSequences(proteinSeqs)
        if len(proteinSeqs) > 1 and proteinSeqs[0] == proteinSeqs[1]:
            loopFound = True

    print("The final alignment is:")
    print(proteinSeqs[0] + "\n")

#run main function
main()