import numpy as np

def configuration(Lr, Lx, Ly, scenario_arr, wheat_fraction):
    """
    Return the spatial arrangement of patches.

    Return the spatial arrangement of patches within the field 
    and the wheat fraction within each patch.

    Parameters
    ----------
    Lr : int 
        number of seasons
    Lx : int 
        number of patch along the x-axis
    Ly : int 
        number of patch along the y-axis
    scenario_rot : str
        rotation scenario (chose: uniform, random, chessboard, alternate, alternate_rank,
                                alternate_strip2, alternate_strip3, alternate_pairs, etc)
    wheat_fraction : float
        wheat fraction within each patch

    Returns
    -------
    array
    """
    arrangement = np.ones((Lr, Lx, Ly))

    ### MIX OF X% OF WHEAT IN EACH PLOT FOR ALL PLOTS
    if scenario_arr == "uniform":
        arrangement = wheat_fraction * arrangement

    ### RANDOM PATTERN OF X% OF WHEAT
    elif scenario_arr == "random":
        id = np.random.permutation(Lx * Ly)
        arrangement = arrangement.reshape([Lr*Lx*Ly])
        arrangement[id[0:int(np.floor(Lx*Ly*(1-wheat_fraction)))]] = 0
        arrangement = arrangement.reshape([Lr, Lx, Ly])

    ### CHESSBOARD PATTTERN WITH 50% OF WHEAT
    elif scenario_arr == "chessboard":
        for i in range(0,Lx):
            for j in range(0,Ly):
                if (i+j)%2 == 0:
                    arrangement[0,i,j] = 0

    ### ALTERNATE ROWS
    elif scenario_arr =="alternate":
        for i in range(0,Lx):
            if i % 2 == 1:
                arrangement[:,i,:] = 0

    ### ALTERNATE RANKS (ALONG Ly)
    elif scenario_arr =="alternate_rank":
        for i in range(0,Ly):
            if i % 2 == 1:
                arrangement[:,:,i] = 0

    ### ALTERNATE STRIPS WITH 2 ROWS PER STRIP
    elif scenario_arr == "alternate_strip2":
        for i in range (0,Lx):
            if i % 4 == 0:
                arrangement[:,i:i+2,:] = 0

    ### ALTERNATE STRIPS WITH 3 ROWS PER STRIP
    elif scenario_arr == "alternate_strip3":
        for i in range (0,Lx):
            if i % 6 == 0:
                arrangement[:, i:i+3, :] = 0

    ### ALTERNATE PAIRS OF PATCH ALONG Lx
    elif scenario_arr == "alternate_pairs":
        for i in range(0,Lx):
            for j in range(0,Ly):
                if i%2 == 0 and j%4 == 0:
                    arrangement[:, i, j:j+2] = 0
                    arrangement[:, i+1, j+2:j+3] = 0

    ### ALTERNATE DOUBLE-PAIRS ALONG Lx
    elif scenario_arr == "alternate_doublepairs":
        for i in range(0,Lx):
            if i % 2 == 1:
                arrangement[:,i,:] = 0
        for i in range(0,Lx):
            for j in range(0,Ly):
                if i%2 == 0 and j%8 == 0:
                    arrangement[:, i, j:j+4] = 0
                    arrangement[:, i+1, j:j+4] = 1

    ### TWO SUB-FIELDS
    elif scenario_arr == "two_subs":
        for i in range(0, Lx):
            if i < (Lx/2):
                arrangement[:, i, :] = 0

    ### FOUR SUB-FIELDS
    elif scenario_arr == "four_subs":
        for i in range(0, Lx):
            for j in range(0, Ly):
                if i < (Lx/2) and j < (Ly/2):
                    arrangement[:, i, j] = 0
                if i >= (Lx/2) and j >= (Ly/2):
                    arrangement[:, i, j] = 0

    ### 4-PATCH SQUARE
    elif scenario_arr == "4_square":
        for i in range(0, Lx):
            for j in range(0, Ly):
                if i % 4 == 0 and j % 4 == 0:
                    arrangement[:, i:i+2, j:j+2] = 0
                    arrangement[:, i+2:i+4, j+2:j+4] = 0

    ### 9-PATCH SQUARES
    elif scenario_arr == "9_square":
        for i in range(0, Lx):
            for j in range(0, Ly):
                if i % 6 == 0 and j % 6 == 0:
                    arrangement[:, i:i+3, j:j+3] = 0
                    arrangement[:, i+3:i+6, j+3:j+6] = 0

    ############ TEST DIFFERENT PROPORTIONS ACCORDING TO ARRANGEMENT ############
    ### 1/5 IN ALTERNATE ROWS
    elif scenario_arr == "1_5_alternate":
        for i in range(0, Lx):
            if i % 5 == 0:
                arrangement[:, i, :] = 0

    ### 2/5 IN ALTERNATE ROWS
    elif scenario_arr == "2_5_alternate":
        for i in range(0, Lx):
            for j in range(0, Lx):
                if i % 5 == 0:
                    arrangement[:, i, :] = 0
                if j % 5 == 2:
                    arrangement[:, j, :] = 0

    ### 3/5 IN ALTERNATE ROWS
    elif scenario_arr == "3_5_alternate":
        for i in range(0, Lx):
            for j in range(0, Lx):
                if i % 5 == 0:
                    arrangement[:, i, :] = 0
                if j % 5 == 2:
                    arrangement[:, j:j+2, :] = 0

    ### 4/5 IN ALTERNATE ROWS
    elif scenario_arr == "4_5_alternate":
        for i in range(0, Lx):
            if i % 5 == 1:
                arrangement[:, i:i+4, :] = 0

    ### 1/5 IN 4-PATCH SQUARE
    elif scenario_arr == "1_5_square":
        for i in range(0, Lx):
            for j in range(0, Ly):
                if i % 10 == 0 and j % 10 == 0:
                    arrangement[:, i:i+2, j:j+2] = 0
                    arrangement[:, i+2:i+4, j+4:j+6] = 0
                    arrangement[:, i+4:i+6, j+8:j+10] = 0
                    arrangement[:, i+6:i+8, j+2:j+4] = 0
                    arrangement[:, i+8:i+10, j+6:j+8] = 0

    ### 2/5 IN 4-PATCH SQUARE
    elif scenario_arr == "2_5_square":
        for i in range(0, Lx):
            for j in range(0, Ly):
                if i % 10 == 0 and j % 10 == 0:
                    arrangement[:, i:i+2, j:j+2] = 0
                    arrangement[:, i:i+2, j+4:j+6] = 0
                    arrangement[:, i+2:i+4, j+2:j+4] = 0
                    arrangement[:, i+2:i+4, j+6:j+8] = 0
                    arrangement[:, i+4:i+6, j+4:j+6] = 0
                    arrangement[:, i+4:i+6, j+8:j+10] = 0
                    arrangement[:, i+6:i+8, j:j+2] = 0
                    arrangement[:, i+6:i+8, j+6:j+8] = 0
                    arrangement[:, i+8:i+10, j+2:j+4] = 0
                    arrangement[:, i+8:i+10, j+8:j+10] = 0

    ### 3/5 IN 4-PATCH SQUARE
    elif scenario_arr == "3_5_square":
        for i in range(0, Lx):
            for j in range(0, Ly):
                if i % 10 == 0 and j % 10 == 0:
                    arrangement[:, i:i+2, j:j+2] = 0
                    arrangement[:, i:i+2, j+4:j+8] = 0
                    arrangement[:, i+2:i+4, j+2:j+4] = 0
                    arrangement[:, i+2:i+4, j+6:j+10] = 0
                    arrangement[:, i+4:i+6, j:j+2] = 0
                    arrangement[:, i+4:i+6, j+4:j+6] = 0
                    arrangement[:, i+4:i+6, j+8:j+10] = 0
                    arrangement[:, i+6:i+8, j:j+4] = 0
                    arrangement[:, i+6:i+8, j+6:j+8] = 0
                    arrangement[:, i+8:i+10, j+2:j+6] = 0
                    arrangement[:, i+8:i+10, j+8:j+10] = 0

    ### 4/5 IN 4-PATCH SQUARE
    elif scenario_arr == "4_5_square":
        for i in range(0, Lx):
            for j in range(0, Ly):
                if i % 10 == 0 and j % 10 == 0:
                    arrangement[:, i:i+2, j:j+8] = 0
                    arrangement[:, i+2:i+4, j+2:j+10] = 0
                    arrangement[:, i+4:i+6, j+4:j+10] = 0
                    arrangement[:, i+4:i+6, j:j+2] = 0
                    arrangement[:, i+6:i+8, j:j+4] = 0
                    arrangement[:, i+6:i+8, j+6:j+10] = 0
                    arrangement[:, i+8:i+10, j:j+6] = 0
                    arrangement[:, i+8:i+10, j+8:j+10] = 0

    return arrangement


