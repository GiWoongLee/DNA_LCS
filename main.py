from consts import *
from sys import argv
import argparse

parser = argparse.ArgumentParser(description='Hairpin Finder by woong. Please Enjoy!')
parser.add_argument('-v', '--verbose', help='increase output verbosity', action='store_true')
parser.add_argument('input_file', help='input file to read a DNA sequence', type=argparse.FileType('r'))
args = parser.parse_args()

# #Main Function
# #read file name from script argument and deliver it to setDNAseq()
def main(DNASeq):
    # Exception - Wrong File Input
    # Exception - LOOP_ITERATION_NUM Too short

    pos = 0
    nextPosAfterHairPin = 0
    #Iterate Through DNA_SEQ

    while(pos< len(DNASeq[0])-MAX_EXP_HAIRPIN_LENGTH):
        for pos in range(nextPosAfterHairPin,len(DNASeq[0])-MAX_EXP_LOOP_LENGTH):
            temp = kmerExist(DNASeq[0],pos)
            if temp['exist'] == True:
                candRegion = getCandidateRegion(DNASeq[0],temp['lKMerStartPos'],temp['lKMerEndPos'],temp['rKMerStartPos'],temp['rKMerEndPos'])
                lLcsInfo = findLCS(candRegion['lsStr'],candRegion['rsStr'][::-1])
                rLcsInfo = findLCS(candRegion['csStr'],candRegion['csStr'][::-1])
                processedLcsInfo1 = processLcsInfo(lLcsInfo[0],lLcsInfo[1],lLcsInfo[2],[],len(candRegion['lsStr']),len(candRegion['rsStr']))
                processedLcsInfo2 = processLcsInfo(rLcsInfo[0], rLcsInfo[1], rLcsInfo[2], [], len(candRegion['csStr']), len(candRegion['csStr']))
                extendInfo = extendKMerToLCS(processedLcsInfo1,candRegion['leftKMer'],processedLcsInfo2)
                isHairPin = finalCheck(processedLcsInfo1,candRegion['leftKMer'],processedLcsInfo2,extendInfo[0],extendInfo[1],extendInfo[2])
                if(isHairPin):
                    # print(candRegion['leftKMer'])
                    # print(candRegion['rightKMer'])
                    # print(temp['lKMerStartPos'])
                    # print(temp['rKMerStartPos'])
                    # print(candRegion['csStr'])
                    # print(candRegion['lsStr'])
                    # print(candRegion['rsStr'])
                    # print(candRegion['rsStr'][::-1])
                    # print(len(candRegion['rsStr'][::-1]))
                    # print(extendInfo)
                    nextPosAfterHairPin = pos + MAX_EXP_HAIRPIN_LENGTH - len(candRegion['lsStr'])
                    break
                else:
                    continue
            else:
                continue

# Input : DNA Sequence
# Output : KMer exist?, Info about KMer
def kmerExist(DNASeq,currPos):
    baseKMer = DNASeq[currPos:currPos + K_MER_LENGTH]

    #Find possible targetKMer position regarding distance
    for distance in range(MAX_EXP_LOOP_LENGTH,MAX_EXP_HAIRPIN_LENGTH+1):
        #Exception
        if(currPos+distance >= len(DNASeq)):
            break

        targetKMer = DNASeq[currPos + distance - K_MER_LENGTH:currPos + distance]
        if baseKMer == targetKMer[::-1]:
            return {'exist':True, 'lKMerStartPos':currPos,'lKMerEndPos':currPos+K_MER_LENGTH-1,'rKMerStartPos':currPos+distance-K_MER_LENGTH,'rKMerEndPos':currPos+distance-1}
        else:
            continue
    return {'exist': False}

# Input : Two DNA Sequences
# Output : LCS Info(LCS length and match info of each letters)
def findLCS(base,target):
    m = len(base)
    n = len(target)
    matchResult = [["" for x in range(m+1)] for y in range(n+1)]
    lcsLength = [[0 for x in range(m+1)] for y in range(n+1)]
    charPair = [[((),) for x in range(m+1)] for y in range(n+1)]

    for i in range(1,m+1):
        lcsLength[i][0] = 0

    for j in range(1,n+1):
        lcsLength[0][j] = 0

    for i in range(1,m+1):
        for j in range(1,n+1):
            if base[i-1] == target[j-1]:
                lcsLength[i][j] = lcsLength[i-1][j-1] + 1
                matchResult[i][j] = "ARROW_CORNER"
            elif lcsLength[i-1][j] >= lcsLength[i][j-1]:
                lcsLength[i][j] = lcsLength[i-1][j]
                matchResult[i][j] = "ARROW_UP"
            else:
                lcsLength[i][j] = lcsLength[i][j-1]
                matchResult[i][j] = "ARROW_LEFT"
            charPair[i][j] = (base[i-1], target[j-1])
    return (matchResult, lcsLength,charPair)


#Input: matchResult, lcsLength and charPair from findLCS && lsStrLcsInfo, rsStrLcsInfo
#Output: lists of tuples (character, matchInfo)
def processLcsInfo(matchResult,lcsLength,charPair,res,i,j):
    if i==0 and j==0:
        return res
    elif i==0:
        for k in range(j,0,-1):
            res.insert(0,("INDEL",('-',charPair[1][k][1])))
        return res
    elif j==0:
        for k in range(i,0,-1):
            res.insert(0,("INDEL",(charPair[k][1][0],'-')))
        return res
    else:
        if(matchResult[i][j]=='ARROW_CORNER'):
            res.insert(0,("MATCH",charPair[i][j]))
            return processLcsInfo(matchResult,lcsLength,charPair,res,i-1,j-1)
        elif(matchResult[i][j]=="ARROW_UP"):
            res.insert(0,("INDEL",(charPair[i][j][0],'-')))
            return processLcsInfo(matchResult, lcsLength, charPair,res, i - 1, j)
        else:
            res.insert(0,("INDEL",('-',charPair[i][j][1])))
            return processLcsInfo(matchResult,lcsLength,charPair,res,i,j-1)

#Input: processedLCSInfo(from lsStr,rsStr'), processedLCSInfo(from msStr,msStr'), leftKMer
#Output: last Match Position within left LCS part and right LCS part
def extendKMerToLCS(lLcsInfo,kMer,rLcsInfo):
    tempPos = len(lLcsInfo)
    tempPos2 = 0
    lStrlist = []
    rStrlist = []

    for pos in range(len(lLcsInfo)-1,-1,-1):
        if lLcsInfo[pos][0] == 'INDEL':
            pass
        else:
            #Extending Direction - left
            if isStopCondition(lLcsInfo,pos,"LEFT"):
                break
            else:
                #update lastMatchPosition
                tempPos= pos
                lStrlist.insert(0,lLcsInfo[pos][1][1])

    for pos in range(len(rLcsInfo)):
        if rLcsInfo[pos][0] == 'INDEL':
            pass
        else:
            # Extending Direction - right
            if isStopCondition(rLcsInfo,pos,"RIGHT"):
                break
            else:
                # update lastMatchPosition
                tempPos2 = pos
                rStrlist.append(rLcsInfo[pos][1][1])

    return ((tempPos,tempPos2), ''.join(lStrlist),''.join(rStrlist))

#Input : Info about candidate hairpin sequence
#Output : is candidate hairpin sequence? + if yes, return lcs and loop string
def finalCheck(lLcsInfo,KMer,rLcsInfo,lastMatchPos,lLcsStr,rLcsStr):
    lcsSeq = lLcsStr + KMer + rLcsStr
    hairpinStr = []
    loopStr = []

    for i in range(lastMatchPos[0],len(lLcsInfo)):
        if(lLcsInfo[i][1][0] == '-'):
            pass
        else:
            hairpinStr.append(lLcsInfo[i][1][0])


    for ch in KMer:
        hairpinStr.append(ch)

    for j in range(lastMatchPos[1]+1):
        if(rLcsInfo[j][1][0]== '-'):
            pass
        else:
            hairpinStr.append(rLcsInfo[j][1][0])

    for i in range(lastMatchPos[1]+1,len(rLcsInfo)-lastMatchPos[1]-1):
        if(rLcsInfo[i][1][0]=='-'):
            pass
        else:
            hairpinStr.append(rLcsInfo[i][1][0])
            loopStr.append(rLcsInfo[i][1][0])

    for j in range(lastMatchPos[1],-1,-1):
        if (rLcsInfo[j][1][1] == '-'):
            pass
        else:
            hairpinStr.append(rLcsInfo[j][1][1])

    for ch in KMer[::-1]:
        hairpinStr.append(ch)

    for i in range(len(lLcsInfo)-1,lastMatchPos[0]-1,-1):
        if(lLcsInfo[i][1][1]== '-'):
            pass
        else:
            hairpinStr.append(lLcsInfo[i][1][1])

    hairpinSeq = ''.join(hairpinStr)
    loopSeq = ''.join(loopStr)

    if(isHairpinConstraint(len(loopSeq),len(hairpinSeq))):
        print(hairpinSeq)
        print(lcsSeq)
        print(loopSeq)
        return True
    else:
        return False


#Helper Functions

#Input : DNA_SEQ, left-k-mer, right-k-mer, length and index
#Output : left-surplus-string, k-mer, center-surplus-string, k-mer', right-surplus-string
def getCandidateRegion(DNASeq,lKMerStartPos,lKMerEndPos,rKMerStartPos,rKMerEndPos):
    csStrLen = rKMerStartPos - lKMerEndPos -1
    lrSurplusStrLen = (MAX_EXP_HAIRPIN_LENGTH - 2 * K_MER_LENGTH - csStrLen)

    #Exception - if there is not enough number of lsStr comparing with lrSurplusStrLen
    if(lKMerStartPos < lrSurplusStrLen//2):
        lsStr = DNASeq[:lKMerStartPos]
    else:
        lsStr = DNASeq[lKMerStartPos-lrSurplusStrLen//2:lKMerStartPos]

    # Exception - if there is not enough number of rsStr comparing with lrSurplusStrLen
    if(len(DNASeq) - rKMerEndPos < lrSurplusStrLen//2):
        rsStr = DNASeq[rKMerEndPos+1:]
    else:
        rsStr = DNASeq[rKMerEndPos+1:rKMerEndPos+lrSurplusStrLen//2+1]

    return {'lsStr':lsStr, 'leftKMer':DNASeq[lKMerStartPos:lKMerEndPos+1],'csStr':DNASeq[lKMerEndPos+1:rKMerStartPos], 'rightKMer':DNASeq[rKMerStartPos:rKMerEndPos+1],'rsStr':rsStr}

#Input : two string subsequences
#Output : boolean(true/false)
def isPalindrome(base,target):
    ###Exception Handling ###
    #base case
    if(len(base)==0 and len(target)==0):
        return True
    #recursive case
    elif(base[0]==target[-1]):
        return isPalindrome(base[1:],target[:-1])
    else:
        return False


#Input : two string related infos(match result, character) / current iterating position
#Output : yes or no
def isStopCondition(iterTarget,pos,extendingDirection):
    if(extendingDirection=="LEFT"):
        #Exception
        if(pos+1>=len(iterTarget)):
            return False

        if(iterTarget[pos+1][0]=="INDEL"):
            return checkIndelNum(iterTarget,pos,extendingDirection)
        else:
            return False
    else:
        # Exception
        if (pos - 1 < 0):
            return False

        if(iterTarget[pos-1][0]=="INDEL"):
            return checkIndelNum(iterTarget, pos, extendingDirection)
        else:
            return False

#Input : iterTarget,pos,extendingDirection
#Output : whether indelNumber is greater than or equal to alpha within the beta position size region
def checkIndelNum(iterTarget,pos,extendingDirection):

    #Exception - Within searchBoundary, KMer Exists
    count = 0
    searchBoundary = (BETA_REGION_LENGTH -1)//2
    for i in range(1,searchBoundary+1):
        #Left Side Indel Checking
        if(pos-i<0):
            pass
        else:
            if(iterTarget[pos-i][0]=="INDEL"):
                count += 1
        #Right Side Indel Checking
        if(pos+i>len(iterTarget)-1):
            pass
        else:
            if (iterTarget[pos + i][0] == "INDEL"):
                count +=1
    return (count >= ALPHA_INDEL_COUNT)


#Input : Info about candidate hairpin sequence
#Output : is candidate hairpin sequence? + if yes, return lcs and loop string
def isHairpinConstraint(loopLength,hairpinSequenceLength):
    return (0 <= loopLength <= MAX_EXP_LOOP_LENGTH) and (MIN_EXP_HAIRPIN_LENGTH <= hairpinSequenceLength <= MAX_EXP_HAIRPIN_LENGTH)

if __name__ == '__main__':
    DNASeq = ""
    with args.input_file as f:
        DNASeq = f.readlines()[1:]
        f.close()

    main(DNASeq)