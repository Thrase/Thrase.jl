# EToV = [1 5 8 2 14 9 16 15 12;
#        8 4 15 6 10 8 1 7 3;
#        3 6 2 11 7 10 9 6 11
#        2 13 6 13 15 15 8 5 2]

# EToF = [20 24 18 19 2 15 4 11 23;
#        18 22 11 14 1 7 8 12 9
#        8 13 7 5 6 3 10 16 21;
#        9 14 5 17 16 1 3 24 19]

# FToB = [0 2 0 1 7 1 0 0 7 2 0 2 2 0 1 0 1 0 0 2 2 1 1 7]

# (FToE, FToLF, EToO, EToS) = connectivityarrays_2d(EToV, EToF)

EToV = [8 4;
        7 6;
        10 3;
        6 2;
        12 5;
        5 7;
        1 9;
        4 11]


EToF = [1 7;
        2 8;
        3 9;
        4 10;
        5 2;
        6 11]

(FToE, FToLF, EToO, EToS) = connectivityarrays_3d(EToV, EToF)
        


EToV = [9 5 15 8;
        7 13 16 9;
        1 2 3 14;
        17 7 4 10;
        10 6 13 3;
        13 12 5 7;
        11 18 7 15;
        12 17 2 13]

EToF = [18 16 6 5;
        17 17 8 2;
        11 9 14 4;
        1 20 9 6;
        2 15 13 3;
        12 10 7 19];

(FToE, FToLF, EToO, EToS) = connectivityarrays_3d(EToV, EToF)

dx1 = [5 1 7 3];
dx2 = [6 2 8 4];
dx3 = [1 2 3 4];
dx4 = [5 6 7 8];
dx5 = [1 2 5 6];
dx6 = [3 4 7 8];

