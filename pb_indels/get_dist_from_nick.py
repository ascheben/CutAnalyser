import sys

def solve(r1, r2):
     # sort the two ranges such that the range with smaller first element
     # is assigned to x and the bigger one is assigned to y
     x, y = sorted((r1, r2))
     #now if x[1] lies between x[0] and y[0](x[1] != y[0] but can be equal to x[0])
     #then the ranges are not overlapping and return the differnce of y[0] and x[1]
     #otherwise return 0 
     if x[0] <= x[1] < y[0] and all( y[0] <= y[1] for y in (r1,r2)):
        return y[0] - x[1]
     return 0

with open(sys.argv[1], 'r') as mastertab:
    for line in mastertab:
        line = line.strip()
        line = line.split("\t")
        indelrange = [int(line[3]),int(line[4])]
        nickrange = [int(line[6]),int(line[7])]
        distance = solve(indelrange,nickrange)
        #if indel downstream, make dist negative
        if int(line[4]) < int(line[6]):
            distance = "-" + str(distance)
        elif distance > 0:
            distance = "+" + str(distance)
        outline = line + [distance]
        print(*outline,sep="\t")
