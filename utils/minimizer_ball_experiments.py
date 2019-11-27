import random
import string
import itertools
import Levenshtein

l = 11
percentage_retain_hashes=0.0005

def randomString(stringLength=10):
    """Generate a random string of fixed length """
    letters = "ACTG"
    return ''.join(random.choice(letters) for i in range(stringLength))

genome=randomString(1000)
#print(genome)

reg_minimizers = set()
space_size = 4**l

for lmer in (''.join(x) for x in itertools.product('ACTG', repeat=l)):
    h = hash(lmer) % space_size
    if h < space_size*percentage_retain_hashes:
        reg_minimizers.add(lmer)

print("kept",len(reg_minimizers),"regular minimizers (%.02f%%)" % (len(reg_minimizers)*100.0/space_size))

def minimizer_spacing(minimizer_set):
    last_seen=0
    positions = []
    for i in range(len(genome)):
        lmer = genome[i:i+l]
        if lmer in minimizer_set:
            positions += [i]
        lmer = genome[i:i+l-1]
        if lmer in minimizer_set:
            positions += [i]
        lmer = genome[i:i+l+1]
        if lmer in minimizer_set:
            positions += [i]

    print(positions)
    if len(positions) == 0: return -1
    mean_spacing = sum([positions[i+1]-positions[i] for i in range(len(positions)-1)])/(1.0*len(positions))
    return mean_spacing

reg_meanspace = minimizer_spacing(reg_minimizers)
print("mean distance between regular minimizers: %.2f" % reg_meanspace)

# ball of diameter 1
def levenshtein_ball(lmer):
    for pos in range(1,l):
        lmer_mutated = list(lmer)
        for c in "ACTG":
            lmer_mutated[pos] = c
            s = ''.join(lmer_mutated)
            if s == lmer: continue
            yield s
    for pos in range(1,l-1):
        lmer_mutated = lmer[:pos] + lmer[pos+1:]
        yield lmer_mutated
    for pos in range(1,l-1):
        for additional in "ACTG":
            lmer_mutated = lmer[:pos] + additional + lmer[pos:]
            yield lmer_mutated



test_lmer=("ACTG"*int(l/4+1))[:l]
for lmer in levenshtein_ball(test_lmer):
    assert(Levenshtein.distance(lmer,test_lmer) == 1 or (len(lmer) == (l+1) and Levenshtein.distance(lmer,test_lmer) == 2))

lev_minimizers = reg_minimizers | set([x for y in reg_minimizers for x in levenshtein_ball(y)])
print("now",len(lev_minimizers),"balled minimizers (%.02f%%)" % (len(lev_minimizers)*100.0/space_size))


lev_meanspace = minimizer_spacing(lev_minimizers)
print("mean distance between balled minimizers: %.2f" % lev_meanspace)

greedy_lev_minims = set()
for reg_min in reg_minimizers:
    can_insert = all((x not in greedy_lev_minims for x in levenshtein_ball(reg_min)))
    if can_insert:
        for x in levenshtein_ball(reg_min):
            greedy_lev_minims.add(x)

print("now",len(greedy_lev_minims),"greedily inserted balled minimizers (%.02f%%)" % (len(greedy_lev_minims)*100.0/space_size))

glev_meanspace = minimizer_spacing(greedy_lev_minims)
print("mean distance between balled minimizers: %.2f" % glev_meanspace)


