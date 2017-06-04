from consts import *
from main import *

import unittest

class TestBasicFunctions(unittest.TestCase):

    def test_getCandidateRegion(self):
        main()

        # print(kmerExist(23))
        # res = for_test_getCandidateRegion(23,32,78,87)
        # print(res)
        # lcsRes1 = findLCS(res['lsStr'],res['rsStr'][::-1])
        # print(lcsRes1)
        # lcsRes2 = findLCS(res['csStr'],res['csStr'][::-1])
        # print(lcsRes2)
        # extendRes1 = processLcsInfo(lcsRes1[0],lcsRes1[1],lcsRes1[2],[],len(res['lsStr']),len(res['rsStr']))
        # print(extendRes1)
        # extendRes2 = processLcsInfo(lcsRes2[0],lcsRes2[1],lcsRes2[2],[],len(res['csStr']),len(res['csStr']))
        # print(extendRes2)
        # final = extendKMerToLCS(extendRes1,res['leftKMer'],extendRes2)
        # print(final)
        # checkHairpin = finalCheck(extendRes1,res['leftKMer'],extendRes2,final[0],final[1],final[2])



if __name__ == '__main__':
    unittest.main()
