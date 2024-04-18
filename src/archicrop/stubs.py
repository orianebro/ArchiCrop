import numpy as np



def leaf_and_internode_appearance(cum_thermal_time, phyllochron, max_nb_of_phytomers, leaves, internodes):
        r"""Triggers the appearance of a new leaf on internode N and of internode N+1, according to the phyllochron.

        

        Parameters
        ----------
        cum_thermal_time : float
            Cumulated thermal time.
        phyllochron : float
            Phyllochron of the species/variety.
        max_nb_of_phytomers : int
            Maximal number of phytomers (phytomer = internode + leaf)
        leaves : list or array
            Sequence of leaves of the growing plant.
        internodes : list or array
            Sequence of internodes of the growing plant.

        Returns
        -------
        leaves : list or array
            Sequence of leaves of the growing plant.
        internodes : list or array
            Sequence of internodes of the growing plant.
        """
        # if it is time for a new leaf and the next internode to grow
        # if cum_thermal_time%phyllochron == 0:
            # new leaf at internode n
            # new internode n+1

        return leaves, internodes



def internode_elongation(internode):
        r"""

        Parameters
        ----------
        internode : list or array
            Sequence of internodes of the growing plant.

        Returns
        -------
        internode : list or array
            Sequence of internodes of the growing plant.
        """

        # in internodes --> growing internodes (bool growing ?)
        # split delta plant height among growing internodes 

        # def plant_height(t):
        #     return 300 / (1 + np.exp(-0.1*(t-100/2)))

        # only one internode growing at a time
        # internodes[-1] += plant_height(t)-plant_height(t-1) # or internodes[i] = plant_height(t)

        return internode



def leaf_elongation(leaf):
        r"""Triggers the appearance of a new leaf on internode N and of internode N+1, according to the phyllochron.

        Parameters
        ----------
        internodes : list or array
            Sequence of internodes of the growing plant.

        Returns
        -------
        internodes : list or array
            Sequence of internodes of the growing plant.
        """

        # in internodes --> growing internodes (bool growing ?)
        # split delta plant height among growing internodes 

        # def plant_height(t):
        #     return 300 / (1 + np.exp(-0.1*(t-100/2)))

        # only one internode growing at a time
        # internodes[-1] += plant_height(t)-plant_height(t-1) # or internodes[i] = plant_height(t)

        return leaf