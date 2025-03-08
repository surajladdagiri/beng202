#  main file to run our code
import Score

fluorescent_seqs = {    
    "mangoI": "GGCACGUACGAAGGGACGGUGCGGAGAGGAGAGUACGUGC",
    "mangoII": "GGCACGUACGAAGGAGAGGAGAGGAAGAGGAGAGUACGUGC",
    "mangoIII": "GGCACGUACGAAGGAAGGAUUGGUAUGUGGUAUAUUCGUACGUGCC",
    "mangoIV": "GGCACGUACCGAGGGAGUGGUGAGGAUGAGGCGAGUACGUGC",
    "beetroot": "GUUAGGCAGAGGUGGGUGGUGUGGAGGAGUAUCUGUC", 
    "spinach": "GACGCGACUGAAUGAAAUGGUGAAGGACGGGUCCAGGUGUGGCUGCUUCGGCAGUGCAGCUUGUUGAGUAGAGUGUGAGCUCCGUAACUAGUCGCGUC", 
    "corn": "GGCGCGAGGAAGGAGGUCUGAGGAGGUCACUGCGCC",
    "chili": "GGCUAGCUGGAGGGGCGCCAGUUCGCUGGUGGUUGGGUGCGGUCGGCUAGCC",
    "pepper": "GGCGCACUGGCGCUGCGCCUUCGGGCGCCAAUCGUAGCGUGUCGGCGCC",
    "mg": "GGUACCCGACUGGCGAGAGCCAGGUAACGAAUGGUACC",
    "tmr3": "GACUCAUUUCCGUUUUCUAUGGGUCUUGGCCUGCUUCGGCAGGAGGUAACAUACACCUGUCCACACCUGCGC"
    }

gene = ""

if not fluorescent_seqs:
    raise RuntimeError("Please pass in some fluorescent sequences!")

if not gene:
    raise RuntimeError("Please pass in a gene sequence!")

best_score = float('-inf')
best_f = ""
for i in range(len(fluorescent_seqs)):
    curr_f = fluorescent_seqs[i]
    sub_str, score = Score.score(gene, curr_f)
    if best_score < score:
        best_score, best_f = score, curr_f
    print(f"For fluorophore sequence {i + 1}:")
    print(f"Optimal substring in gene: {sub_str}\tScore: {score}")

print(f"Best fluorophore: {best_f}\nScore: {best_score}")
