
f = open("C:/Users/Gunnar/Documents/Ph.D/Learning/DMRG/Tryout code/TEBD/tester.txt","r")
inputlist = []
dic = {}
for line in f:
    if line.startswith('#'):
        continue
    inputlist.append(line.rstrip('\n').split(" = "))
    x = inputlist[-1][0]
    y = inputlist[-1][1]
    try:
        d = int(y)
        
    except ValueError:
        try:
            d = float(y)
        except ValueError:
            d = y
    finally:
        dic[x] = d
f.close()
