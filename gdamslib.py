#!/usr/bin/python
# -*- coding: utf-8 -*-
""" Design of small gabion dams for soil conservation

Copyright (c) 2014-2018 Eduardo Jiménez <ecoslacker@irriapps.com>

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import sys
import csv
import math
import codecs
import pprint
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
# matplotlib.use('Qt4Agg')

# Define constants
LENGTH = 0
WIDTH = 1
HEIGHT = 2
COORDX = 3
COORDY = 4
DATA_COLS = 2
MIN_DATA = 3

# CLASSES START


class DataFile:
    """ Management for stream's cross section data files. """

    def __init__(self, filename):
        self.filename = filename
        self.data = []

    def read_data(self):
        with codecs.open(self.filename, 'r', 'utf-8') as f:
            reader = csv.reader(f)
            for row in reader:
                self.data.append(row)

    def check_data(self):

        status = True

        if len(self.data) < MIN_DATA:
            print("The data should has 3 values at minimum")
            return False

        # Check the size of the data array (two columns)
        for row in self.data:
            print(row)
            if len(row) != DATA_COLS:
                print("Length of data arrays does not match")
                return False
            try:
                # Try to convert data to float to catch no numeric values
                float(row[0])
                float(row[1])
            except ValueError:
                print("Could not convert string to float: {0}".format(row))
                return False
        print("Data array length and values are OK!")
        return status

    def show_data(self):
        pprint.pprint(self.data)


class ChannelSection:
    """ Channel cross section

    Put data in useful lists and draw the cross section.
    """

    def __init__(self, data, layer_height):
        self.data = data
        self.data_x = []
        self.data_y = []
        self.layers_height = layer_height

        # Initialize some variables
        self.levels = -1
        self.streambed = -1
        self.maxdepth = -1
        self.channel_width = -1
        self.depths = []
        self.widths = []

    def is_sorted(self, l):
        return all(l[i] <= l[i+1] for i in xrange(len(l)-1))

    def get_xy(self):
        """ Create two list, one for X data and other for Y data. """

        for row in self.data:
            self.data_x.append(float(row[0]))
            self.data_y.append(abs(float(row[1])))

        if not self.is_sorted(self.data_x):
            print("Data is not sorted")
            return False
        else:
            return True

    def show_xy(self):
        """ Show the data lists of X and Y data. """

        pprint.pprint(self.data_x)
        pprint.pprint(self.data_y)

    def draw_channel(self, path=''):
        """ Draw the channel's cross section. """

        plt.plot(self.data_x, self.data_y)
        plt.xlabel(u'Width (m)')
        plt.ylabel(u'Depth (m)')
        plt.title(u'Channel cross section')
        plt.grid(True)
        plt.gca().invert_yaxis()

        # Save fig if there is a path
        if path == '':
            plt.show()
        else:
            plt.savefig(path)
            plt.show()

    def get_levels(self, h):
        """ Get the levels of the channel using the layer's height. """

        # An extra level must be considered for the base,
        # NOTICE: "self.levels" is NOT the number of layers !!
        # it can be two layers in the same level, e.g. both spillway sides

        self.levels = int(self.maxdepth / h) + 1
        return self.levels

    def get_channel_width(self):
        self.channel_width = self.data_x[-1] - self.data_x[0]
        return self.channel_width

    def get_maxdepth(self):
        """ Get the maximum value of depth of channel's cross section. """

        self.maxdepth = max(self.data_y)
        return self.maxdepth

    def get_streambed(self):
        """ Get the center of streambed. """

        if len(self.widths[0]) == 2:
            self.streambed = self.widths[0][0][0]
        else:
            # Assuming two pairs of coordinates
            diff = (abs(self.widths[0][1][0] - self.widths[0][0][0])) / 2.0

            # The stream bed should be in the center
            self.streambed = self.widths[0][0][0] + diff
        return self.streambed

    @staticmethod
    def find_point(x1, x2, y1, y2, d):
        """ Find a point using similar triangles

        Given x1, y1 and x2,y2 finds the "x" coordinate for a given "y"
        between y1 and y2, using similar triangles.
        """

        # Processing large triangle
        base = abs(x2 - x1)
        height = abs(y2 - y1)

        # Put false coordinates to initialize
        x = -1
        y = -1

        if y2 == d:
            x = x2
            y = y2
        elif y2 < d:
            # Second point is deeper than current depth
            if y1 > d:
                # Now processing the small triangle
                h = abs(y2 - d)
                b = (base / height) * h
                x = x2 - b
                y = y2 + h
        elif y2 > d:
            # First point is deeper than current depth
            if y1 < d:
                # Now processing the small triangle
                h = abs(y1 - d)
                b = (base / height) * h
                x = x1 + b
                y = y1 + h
        return x, y

    def get_widths(self):
        """ Obtains widths of the channel

        Return the width of the cross section channel for each
        depth value given.
        """

        self.get_maxdepth()
        self.get_levels(self.layers_height)
        self.get_channel_width()

        # Start-end coordinates for each layer
        self.widths = [0] * self.levels

        for i in range(self.levels):
            current_depth = self.maxdepth - (self.layers_height * i)
            self.depths.append(current_depth)
            # print("   Current depth %d: %.3f" % (i, current_depth))
            p = []  # Start-end coordinates for this layer
            for n in range(len(self.data_x)):
                # WARNING: Check if using n-1 is convenient, it could be bad!
                coordx, coordy = self.find_point(self.data_x[n-1],
                                                 self.data_x[n],
                                                 self.data_y[n-1],
                                                 self.data_y[n],
                                                 current_depth)
                if (coordx == -1) and (coordy == -1):
                    # False coordinates detected, ignore them
                    pass
                else:
                    pair = [coordx, coordy]
                    p.append(pair)
            self.widths[i] = p
            try:
                # Obtaining the width for each depth, using the extreme
                # pair of coordinates
                width = self.widths[i][-1][0] - self.widths[i][0][0]
            except TypeError:
                print("TypeError: Assuming one single point, width = 0")
                width = 0
            except IndexError:
                # Very common, there is only one coordinate pair
                print("IndexError: Assuming one single point, width = 0")
                width = 0
            self.widths[i].append(width)
        # print("Widths array is: %r" % widths)

        # Now widths exist and get streambed
        self.get_streambed()
        return self.widths

    def show_widths(self):
        """ Show widths

        Display the widths of channel cross section at each depth
        using layer's height.
        """

        print(" - Width table (in meters):")
        print("  ---------------------------")
        print("   Level\tDepth\tWidth")
        print("  ---------------------------")
        for i in range(len(self.widths)):
            print("   {0}\t{1:0.3f}\t{2:0.3f}".format(i, self.depths[i],
                                                      self.widths[i][-1]))
        print("  ---------------------------")
        # pprint.pprint(self.widths)

    def channel_report(self):
        # text = "\nCHANNEL CROSS SECTION\n" +\
        #        " Layers height:  {0:0.2f} m\n".format(self.layers_height) +\
        #        " Channel levels: {0}\n".format(self.levels) +\
        #        " Width:          {0:0.2f} m\n".format(self.channel_width) +\
        #        " Maximum depth:  {0:0.2f} m\n".format(self.maxdepth) +\
        #        " Streambed at:   {0:0.2f} m\n".format(self.streambed)
        text = "\nCHANNEL CROSS SECTION\n" + \
               " %-20s %0.2f m\n" % ("Layers height:", self.layers_height) + \
               " %-20s %d\n" % ("Channel levels:", self.levels) + \
               " %-20s %0.2f m\n" % ("Width:", self.channel_width) + \
               " %-20s %0.2f m\n" % ("Maximum depth:", self.maxdepth) + \
               " %-20s %0.2f m\n" % ("Stream bed at:", self.streambed)
        return text


class RationalFormula:
    """ Rational formula

    Use the rational formula to obtain the maximum runoff of the watershed.
    """

    def __init__(self, c=0.0, i=0.0, a=0.0):
        self.coefficient = float(c)   # Runoff coefficient [dimensionless]
        self.rainfall = float(i)      # Rainfall intensity [mm h⁻¹]
        self.area = float(a)          # Area of the watershed [ha]
        self.runoff = None            # Watershed runoff [m³ s⁻¹]

    def get_runoff(self):
        self.runoff = (self.coefficient * self.rainfall * self.area) / 360.0
        return self.runoff


class HighWaterMarks(ChannelSection):
    """ High water marks

    Use the high water marks method to get the maximum runoff of the stream.
    """

    def __init__(self, water_height):
        self.water_height = water_height


class Stream:
    """ Stream

    Stream characteristics, like length, slope, etc.
    """

    def __init__(self, rise, run):
        self.rise = float(rise)
        self.run = float(run)
        self.slope = 0

    def get_slope(self):
        self.slope = self.rise / self.run
        return self.slope


class RunoffCoefficient:
    """ Runoff coefficient

    Get the runoff coefficient of the watershed with land use, soil type,
    and precipitation data. This method is from Official Mexican Norm:
    NOM-011-CONAGUA-2015
    """

    def __init__(self, landuse, soiltype, precipitation):

        self.land_use = int(landuse)
        self.soil_type = int(soiltype)
        self.precipitation = float(precipitation)

        # From this table we can obtain a factor named "k", that is
        # needed to runoff coefficient. Land use and soil type are needed.
        self.table = [[0.26, 0.28, 0.30], [0.24, 0.27, 0.30],
                      [0.24, 0.27, 0.30], [0.24, 0.27, 0.30],
                      [0.14, 0.20, 0.28], [0.20, 0.24, 0.30],
                      [0.24, 0.28, 0.30], [0.07, 0.16, 0.24],
                      [0.12, 0.22, 0.26], [0.17, 0.26, 0.28],
                      [0.22, 0.28, 0.30], [0.26, 0.29, 0.32],
                      [0.27, 0.30, 0.33], [0.18, 0.24, 0.30]]
        self.k = -1
        self.lu = []
        self.st = []
        self.coefficient = -1

    def show_land_use(self):
        self.lu = [u"0.  Barbecho, áreas incultas y desnudas",
                   u"1.  Cultivos en hilera",
                   u"2.  Legumbres o rotación de pradera",
                   u"3.  Granos pequeños",
                   u"4.  Pastizal (75%)",
                   u"5.  Pastizal (50-75%)",
                   u"6.  Pastizal (50%)",
                   u"7.  Bosque (75%)",
                   u"8.  Bosque (50-75%)",
                   u"9.  Bosque (25-50%)",
                   u"10. Bosque (25%)",
                   u"11. Zona urbana",
                   u"12. Caminos",
                   u"13. Pradera permanente"]
        self.lu = ["0.  Bare ground",
                   "1.  Crops",
                   "2.  Vegetables or temporal meadow",
                   "3.  Grain crops",
                   "4.  Grassland (75%)",
                   "5.  Grassland (50-75%)",
                   "6.  Grassland (50%)",
                   "7.  Forest (75%)",
                   "8.  Forest (50-75%)",
                   "9.  Forest (25-50%)",
                   "10. Forest (25%)",
                   "11. Urban",
                   "12. Pathways",
                   "13. Meadow"]
        pprint.pprint(self.lu)

    def show_soil_type(self):

        # self.st = [u"0 = Suelo permeable",
        #            u"1 = Suelo de permeabilidad media",
        #            u"2 = Suelo casi impermeable"]

        self.st = ["Pervious soil", "Semi-pervious soil", "Impervious soil"]
        pprint.pprint(self.st)

    def show_table(self):
        pprint.pprint(self.table)

    def get_coefficient(self):

        # Selecting the k factor for runoff coefficient
        self.k = self.table[self.land_use][self.soil_type]

        if self.k <= 0.15:
            self.coefficient = self.k * (self.precipitation-250.0) / 2000.0
        elif self.k > 0.15:
            self.coefficient = self.k * (self.precipitation-250.0) / 2000.0 \
                               + (self.k - 0.15) / 1.5
        return self.coefficient

    def calc_runoff_coefficient(k, precipitation):
        coefficient = -1
        if k <= 0.15:
            coefficient = k * (precipitation-250.0) / 2000.0
        elif k > 0.15:
            coefficient = k * (precipitation-250.0) / 2000.0 + (k-0.15) / 1.5
        return coefficient
    # The calculation of the runoff coefficient can be static
    calc_runoff_coefficient = staticmethod(calc_runoff_coefficient)


class ConcentrationTime:
    """ Concentration time

    Concentration time of the watershed, is the time it takes the rain from
    the farthest point on watershed to reach the outlet.
    """

    def __init__(self, area, length, rise):
        self.area = float(area)
        self.length = float(length)
        self.rise = float(rise)
        self.concentration_time = -1

    def get_concentration_time(self):
        if self.area >= 200:
            # print(" - (Using Kirpich equation for concentration time)")
            # Use Kirpich equation
            self.concentration_time = ((0.000325 * (self.length ** 0.77) /
                                        ((self.rise/self.length)**0.385))*60)
        elif self.area < 200:
            # Use Giandotti equation
            # Remember: 100 ha = 1 km²
            # print(" - (Using Giandotti equation for concentration time)")
            # print("Area: %.3f ha = %.3f km²" % (self.area, self.area/100))
            # print("Length of main stream: %.3f " % self.length)
            t1 = 4 * math.sqrt(self.area / 100.0)
            t2 = 1.5 * (self.length / 1000.0)
            t3 = math.sqrt((self.rise / self.length) * (self.length / 1000.0))
            self.concentration_time = (((t1 + t2) / (25.3 * t3)) * 60.0)
        return self.concentration_time


class Weir:
    """ Weir

    Hydraulic attributes and operations to dimensioning a broad crested
    rectangular weir, that will be used to dimension the spillway of the dam.
    """

    def __init__(self, q=0, c=1.45, l=0, h=0, n='h', lh=0.5):
        self.flowrate = q     # flow rate
        self.weir_coefficient = c  # weir coefficient for structure
        self.weir_length = l       # width of the crest
        self.weir_head = h         # height of water head over the crest
        self.notch_shape = n  # notch shape
        self.layers_height = lh
        self.weir_height = None
        self.weir_levels = None
        self.weir_layers = None
        self.get_notch()
        # Elements of the weir in the dimensions of the dam
        self.spillway_left_limit = None
        self.spillway_right_limit = None
        self.spillway_left = None
        self.spillway_right = None

        self.weir_notch = -1

    def get_notch(self):
        """ Get exponent: 'h' for horizontal weir or 'v'-notch. """

        if self.notch_shape == 'h':
            self.weir_notch = 3. / 2.
        elif self.notch_shape == 'v':
            self.weir_notch = 5. / 2.

    def get_flowrate(self):
        """ Get the value of flow rate passing through weir. """

        self.get_notch()

        self.flowrate = self.weir_coefficient * self.weir_length * \
            pow(self.weir_head, self.weir_notch)
        return self.flowrate

    def get_weir_length(self):
        """ Get the weir length given head of water and flow rate. """
        self.get_notch()
        self.weir_length = pow(self.flowrate / (self.weir_coefficient *
                               self.weir_head), 1.0 / self.weir_notch)
        return self.weir_length

    def get_weir_head(self):
        """ Get the height of water head over the crest. """

        self.get_notch()

        self.weir_head = pow(self.flowrate / (self.weir_coefficient *
                             self.weir_length), 1.0 / self.weir_notch)
        return self.weir_head

    def get_weir_height(self):
        """ Get the actual weir height according the layer's height used. """

        self.get_notch()
        self.get_weir_head()

        # Weir height
        self.weir_height = 0
        self.weir_levels = 0  # Levels of weir (has 2 layers per level)
        while self.weir_height <= self.weir_head:
            self.weir_levels += 1
            self.weir_height += self.layers_height
        self.weir_layers = self.weir_levels
        return self.weir_height

    def design_weir(self):
        """ Design weir

        Creates the dimensions of the weir, this is actually the same of
        getting the weir height.
        """

        return self.get_weir_height()

    def weir_report(self):
        text = "\nWEIR DESIGN:\n" +\
               " Flow rate:     {0:0.2f} m^3 s^-1\n".format(self.flowrate) +\
               " Head of water: {0:0.2f} m\n".format(self.weir_head) +\
               " Length:        {0:0.2f} m\n".format(self.weir_length) +\
               " Height:        {0:0.2f} m\n".format(self.weir_height) +\
               " Levels:        {0:0.2f}\n".format(self.weir_levels) +\
               " Coefficient:   {0:0.2f}\n".format(self.weir_coefficient) +\
               " Notch shape:   {0}\n".format(self.notch_shape) +\
               " Exponent (n):  {0:0.2f}\n".format(self.weir_notch)

        return text


class Stability(Weir):
    """ Stability

    A simple method for stability analysis as proposed by Oropeza-Mota &
    Lopez-Martinez (2009)

    Stability analysis for small gabion dams for soil conservation
    This assumes the following:
    - Dam is a single body
    - Hydrostatic force is generated by the water and sediments retained
      by the dam body
    - All the operations are performed over an unitary section of the dam
      which exclude the weir layers
    """

    def __init__(self, flowrate, weir_coef, weir_length,
                 water_height, notch, hlayer):
        self.weir = Weir.__init__(self, flowrate, weir_coef, weir_length,
                                  water_height, notch, hlayer)

        # Symbols of dam elements, according to original nomenclature
        self.dam_base = None           # (B)
        self.dam_crest = None          # (b) Length of the top weir layer
        self.dam_volume = None         # (V) Volume of dam (unitary section)
        self.volume_sum = None         # ( ) Sum of (hi * bi) * bi/2, req. Zp
        self.effective_height = None   # (H) Exclude weir and the base layers
        self.dam_weight = None         # (P) Weight of unitary section
        self.dam_lever = None          # (Zp)
        self.water_weight = None       # (q)
        self.water_lever = None        # (Xq)
        self.hydrostatic = None        # (E)
        self.hydrostatic_lever = None  # (XE)
        self.front_lever = None        # (Xp = B -Zp)
        self.sliding = None
        self.turn = None
        self.central = None
        self.total_displacement = None
        self.middle_third = None
        self.dimensions = None  # Dam dimensions table

        self.w = 1.2      # w = specific gravity of water with sediment (t/m^3)
        self.delta = 2.4  # delta = specific gravity of stone (t/m^3)
        self.f = 0.75     # Friction factor for stone
        self.SF = 1       # Security factor

    def analyze(self):
        """ Analyze

        Performs stability analysis
        """
        print("Starting stability analysis of the dam")

        # (B) Base length
        self.dam_base = self.dimensions[0][LENGTH]

        # TODO: Check this!
        # If the weir have more than one layer, is because is required by the
        # water head, so all the layers of the weir should have the same length
        # Otherwise this will be incorrect!
        self.dam_crest = self.dimensions[-1][LENGTH]  # (b) Top of dam length

        self.dam_volume = 0  # (V)
        self.volume_sum = 0  # Sum of volumes and positions of all layers
        self.effective_height = 0  # (H)
        for i in range(len(self.dimensions) - self.weir_layers):
            # Volume of the current layer
            layer_vol = self.dimensions[i][LENGTH] * self.dimensions[i][HEIGHT]

            # Increase the total volume
            self.dam_volume += layer_vol

            # Volume of current layer times it's gravity center position
            self.volume_sum += layer_vol * (self.dimensions[i][LENGTH] / 2.0)

            # Effective height of the dam (H)
            if i is not 0:
                self.effective_height += self.dimensions[i][HEIGHT]

        # Calculate the value of all the elements
        self.get_water_weight()       # (q)
        self.get_water_lever()        # (Xq)
        self.get_dam_weight()         # (P)
        self.get_dam_lever()          # (Zp)
        self.get_hydrostatic()        # (E)
        self.get_hydrostatic_lever()  # (XE)

        # Check conditions
        if not self.check_sliding():
            print(" - Sliding condition not passed!")
            return False

        if not self.check_overturning():
            print(" - Overturning condition not passed!")
            return False

        if not self.check_central():
            print(" - Central condition not passed!")
            return False

        # Verify the central condition
        if not self.check_middle_third():
            print(" - Resultant force is not in middle third!")
            return False

        print(" - All the stability conditions are OK!")
        print("Finishing stability analysis... done!")
        return True

    def get_water_weight(self):
        """ Weight of the water volume over the weir """
        self.water_weight = self.weir_head * self.dam_crest * self.w
        return self.water_weight

    def get_water_lever(self):
        self.water_lever = self.dam_crest / 2.0
        return self.water_lever

    def get_dam_weight(self):
        # Apparent specific gravity = delta - w
        self.dam_weight = self.dam_volume * (self.delta - self.w)
        return self.dam_weight

    def get_dam_lever(self):
        # Please read documentation!
        self.dam_lever = self.volume_sum / self.dam_volume
        return self.dam_lever

    def get_hydrostatic(self):
        # Considering, h = effective_height / 2 and wet area = H * 1
        self.hydrostatic = 1.5 * self.effective_height * self.w
        return self.hydrostatic

    def get_hydrostatic_lever(self):
        self.hydrostatic_lever = self.effective_height / 3.0
        return self.hydrostatic_lever

    def check_sliding(self):
        self.sliding = ((self.water_weight + self.dam_weight) * self.f) / \
            self.hydrostatic
        # The value should be greater or equal than the security factor
        return self.sliding >= self.SF

    def check_overturning(self):
        self.front_lever = self.dam_base - self.dam_lever
        self.turn = (self.dam_weight * self.front_lever) / (self.hydrostatic *
                    self.hydrostatic_lever)
        # The value should be greater or equal than the security factor
        return self.turn >= self.SF

    # TODO: Check this equation! maybe be incorrect or innecessary to use it
    def check_central(self):
        """ Central condition

        The resulting vector of all the forces over the dam, should
        be at 2/3 of the base length """
        self.central = self.water_weight * self.water_lever + \
            self.dam_weight * self.dam_lever + \
            self.hydrostatic * self.hydrostatic_lever
        third = (2. / 3.) * (self.water_weight + self.dam_weight) * \
            self.dam_base
        return self.central <= third

    def check_middle_third(self):
        """ Check middle third

        Resultant force should be located in the middle third of the dam base.
        """
        # NOTE: Check nomenclature for this!
        # Using similar triangles for the resultant force (R) with P and E, let
        # alpha be the angle of R:
        #
        # tan(alpha) =  E / P
        # Also:
        # tan(alpha) = Z' / XE
        #
        # where: Z' is the horizontal displacement of P when E is applied
        #        XE is the lever of the hydrostatic pressure
        #
        # Then:
        # Z' = tan(alpha) * XE

        tan_alpha = self.hydrostatic / self.dam_weight
        # Horizontal displacement (Z')
        horizontal_displacement = tan_alpha * self.hydrostatic_lever

        # Total displacement = Zp + Z', must be less than middle third
        self.total_displacement = self.dam_lever + horizontal_displacement
        self.middle_third = 2 * (self.dam_base / 3.0)
        return self.middle_third > self.total_displacement

    def stability_report(self):
        """ Stability report

        Creates and show a report of the data, values and calculations
        performed during the stability analysis of the dam.
        """

        text = "\nSTABILITY ANALYSIS\n" +\
            "Dimensions:\n" +\
            " Dam effective height (H):      {0:0.2f} m\n".format(self.effective_height) +\
            " Dam base length (B):           {0:0.2f} m\n".format(self.dam_base) +\
            " Crest of dam (b):              {0:0.2f} m\n".format(self.dam_crest) +\
            " Water height (h'):             {0:0.2f} m\n".format(self.weir_head) +\
            " Volume of unitary section (V): {0:0.2f} m^3\n".format(self.dam_volume) +\
            " Wet area (S):                  {0:0.2f} m^2\n".format(self.effective_height) +\
            "\n Forces:\n" +\
            " Water weight (q):            {0:0.2f} t\n".format(self.water_weight) +\
            " Dam weight (P):              {0:0.2f} t\n".format(self.dam_weight) +\
            " Hydrostatic pressure (E):    {0:0.2f} t\n".format(self.hydrostatic) +\
            " Water weight lever (Xq):     {0:0.2f} m\n".format(self.water_lever) +\
            " Dam weight lever (Zp):       {0:0.2f} m\n".format(self.dam_lever) +\
            " Hydrostatic p. lever (XE):   {0:0.2f} m\n".format(self.hydrostatic_lever) +\
            " Dam weight front lever (Xp): {0:0.2f} m\n".format(self.front_lever) +\
            "\n Parameters:\n" + \
            " Specific gravity water-sediment: {0:0.2f} t m^-3\n".format(self.w) + \
            " Specific gravity of rock:        {0:0.2f} t m^-3\n".format(self.delta) +\
            " Apparent dam specific gravity:   {0:0.2f} t m^-3\n".format(self.delta - self.w) +\
            " Friction factor for rock:        {0:0.2f}\n".format(self.f) +\
            "\n Stability conditions:\n" +\
            " Sliding value:            {0:0.2f}\n".format(self.sliding) +\
            " Overturning value:        {0:0.2f}\n".format(self.turn) +\
            " Middle third of the base: {0:0.2f} m\n".format(self.middle_third) +\
            " Resulting force position: {0:0.2f} m\n".format(self.total_displacement)
        return text


class GabionDam(ChannelSection, Stability):
    """ Body of the dam. """

    def __init__(self, section_data, layer_height, flowrate, weir_coef,
                 weir_length, water_height, notch):
        ChannelSection.__init__(self, section_data, layer_height)
        Stability.__init__(self, flowrate, weir_coef, weir_length,
                           water_height, notch, layer_height)

        # Table of dimensions of the dam
        self.abutment = None

    def design(self):
        """ Design

        Generates automatic dimensions for each dam's layer using widths
        list, layer's height, abutment and weir length.
        """

        # Dimension the weir that will be used to create spillway
        self.design_weir()

        # Use the abutment as same as layer's height
        self.abutment = self.layers_height
        print(" - Abutment: {0}".format(self.abutment))

        print("Dimension layers per level:")

        # Automatic generated dimensions list: Length, Width, Height, X, Y
        auto_dimensions = []
        for i in range(self.levels):
            print("LEVEL {0}".format(i))

            # The dam is going to have a 'square' form, it is an equal
            # length and height, so get their maximum dimension
            max_dim = self.levels * self.layers_height
            # Current depth of each level
            d = self.maxdepth - (self.layers_height * i)
            # Set the length value of each layer
            l = max_dim - (self.layers_height * i)
            # Set the width of each layer
            w = self.widths[i][-1] + (self.abutment * 2)
            # Set a constant value for each layer
            h = self.layers_height
            # Set the (x,y) coordinates of left bottom vertex (in front view)
            x = self.widths[i][0][0] - self.abutment
            y = self.widths[i][0][1] + self.layers_height
            # y = d + self.layers_height

            print(" - Current depth: {0:0.2f}".format(d))
            print(" - Length (m):    {0:0.2f}".format(l))
            print(" - Width (m):     {0:0.2f}".format(w))
            print(" - Height (m):    {0:0.2f}".format(h))
            print(" - Coordinate X:  {0:0.2f}".format(x))
            print(" - Coordinate Y:  {0:0.2f}".format(y))

            # Finally append the automatic generated dimensions
            auto_dimensions.append([l, w, h, x, y])

        self.dimensions = auto_dimensions

        # Modify the dam dimensions to create the spillway
        if self.create_spillway():
            print("Weir created successfully!")

        # Adjust all layers to spillway
        self.adjust_dimensions()

        return self.dimensions

    def check_spillway(self):
        """ Check spillway

        Checks if the spillway dimensions are inside the permissible limits
        of the channel and layers dimensions.
        """

        print("Checking spillway...")

        # Get the possible limits of the spillway, this get the layer that
        # the spillway bottom can reach.
        spillway_bottom = self.levels - (self.weir_levels + 1)
        print(" - Dam levels:  {0}".format(self.levels))
        print(" - Weir levels: {0}".format(self.weir_levels))

        print(" - Spillway bottom layer: {0}".format(spillway_bottom))
        # pprint.pprint(self.widths)

        self.spillway_left_limit = self.widths[spillway_bottom][0][0]
        self.spillway_right_limit = self.widths[spillway_bottom][-2][0]

        print(" - Spillway left limit:  {0:0.2f}".format(
                self.spillway_left_limit))
        print(" - Spillway right limit: {0:0.2f}".format(
                self.spillway_right_limit))

        # The spillway should be located at the center of the streambed,
        # or at least cover it.
        self.spillway_left = self.streambed - (self.weir_length / 2.0)
        self.spillway_right = self.spillway_left + self.weir_length

        # Check spillway limits
        print("   Left side: {0}>{1:0.2f}, right side: {2}<{3:0.2f}".format(
                self.spillway_left, self.spillway_left_limit,
                self.spillway_right, self.spillway_right_limit))
        if (self.spillway_left < self.spillway_left_limit or
                self.spillway_right > self.spillway_right_limit):
            print("   Error: Spillway dimensions exceed permissible limits!")
            return False
        else:
            print("Checking spillway limits... OK!")

        # Spillway is inside the permissible limits
        return True

    def create_spillway(self):
        """ Create spillway

        Creates a spillway with the dimensions of the (broad crested)
        rectangular weir, in the upper layers of the dam.
        """

        if not self.check_spillway():
            print("Something went wrong with the spillway!")
            return False

        # WARNING: This will modify the dimensions list
        right_layers = []  # New layers for right side
        start_lvl = self.levels - self.weir_levels
        end_lvl = self.levels
        for level in range(start_lvl, end_lvl):
            print("Creating spillway in level: {0}".format(level))

            # Get width of right spillway layer, before shrink current layer
            width = self.dimensions[level][COORDX] + \
                self.dimensions[level][WIDTH] - self.spillway_right

            # Shrink layer's width to create the spillway's left side
            self.dimensions[level][WIDTH] = self.spillway_left - \
                self.dimensions[level][COORDX]

            # Create the other dimensions of the right side layer
            length = self.dimensions[level][LENGTH]
            height = self.dimensions[level][HEIGHT]
            x = self.spillway_right
            y = self.dimensions[level][COORDY]

            right_layer = []  # Right side spillway layer
            right_layer.append(length)
            right_layer.append(width)
            right_layer.append(height)
            right_layer.append(x)
            right_layer.append(y)

            # Save the current layer
            right_layers.append(right_layer)

        # Add the new layers to the dimensions list
        self.weir_layers = self.weir_levels * 2
        self.dimensions += right_layers

        return True

    def adjust_dimensions(self):
        """ Adjust dimensions

        Modify the layers of the dam as needed to cover the spillway
        """

        for i in range(len(self.dimensions) - self.weir_layers):

            right = self.dimensions[i][COORDX] + self.dimensions[i][WIDTH]

            # Make sure the left side of spillway is covered by all the layers
            if self.dimensions[i][COORDX] > self.spillway_left:
                # Modify layers that not cover the left side of the spillway
                print(" - Spillway not cover by layer {0}: {1} > {2}".format(i,
                      self.dimensions[i][COORDX],
                      self.spillway_left))

                new_width = self.dimensions[i][WIDTH] + \
                    (self.dimensions[i][COORDX] - self.spillway_left)

                print(" -- Changing width from: {0:0.2f} to {1:0.2f}".format(
                        self.dimensions[i][WIDTH], new_width))
                print(" -- Coordinate X from: {0:0.2f} to {1:0.2f}".format(
                        self.dimensions[i][COORDX], self.spillway_left))

                self.dimensions[i][WIDTH] = new_width
                self.dimensions[i][COORDX] = self.spillway_left
            # Make sure the right side of spillway is covered by all the layers
            if right < self.spillway_right:
                # Modify layers that not cover the right side of the spillway
                print(" - Spillway not cover by layer {0}: {1} > {2}".format(
                        i, right, self.spillway_right))

                new_width = self.dimensions[i][WIDTH] + \
                    (self.spillway_right - right)

                print(" -- Changing width from: {0:0.2f} to {1:0.2f}".format(
                        self.dimensions[i][WIDTH], new_width))
                self.dimensions[i][WIDTH] = new_width
        return self.dimensions

    def show_dimensions(self):
        """ Show the dimensions of the dam. """

        # pprint.pprint(self.dimensions)
        print(" - Dimensions table (in meters):")
        print("  -----------------------------------------------------------")
        print("   Layer\t\tLength\tWidth\tHeight\tX\tY")
        print("  -----------------------------------------------------------")
        for i in range(len(self.dimensions)):
            print("   {0}\t{1:0.2f}\t{2:0.2f}\t{3:0.2f}\t{4:0.2f}\t{5:0.2f}".format(
                    i, self.dimensions[i][LENGTH], self.dimensions[i][WIDTH],
                    self.dimensions[i][HEIGHT], self.dimensions[i][COORDX],
                    self.dimensions[i][COORDY]))
        print("  -----------------------------------------------------------")
        return True

    def report_dimensions(self, delimiter="\t"):
        headers = ["Layer", "Length", "Width", "Height", "Coordinate X",
                   "CoordinateY"]
        text = "\nDIMENSIONS OF THE DAM:\n" + delimiter.join(headers) + "\n"
        for i in range(len(self.dimensions)):
            row = []
            for value in self.dimensions[i]:
                row.append("{0:0.2f}".format(value))
            text += "{0}".format(i) + delimiter + delimiter.join(row) + "\n"
        return text

    def front_view_in(self, x=None, y=None, filename='', fig=None, plt=None):

        """ Draws the front view of the dam. """

        if fig is None:
            fig = plt.figure()

        for i in range(len(self.dimensions)):
            left = self.dimensions[i][COORDX]
            right = self.dimensions[i][COORDX] + self.dimensions[i][WIDTH]
            # Because the inverted axis bottom and top should be as follow
            top = self.dimensions[i][COORDY]
            bottom = self.dimensions[i][COORDY] - self.dimensions[i][HEIGHT]

            verts = [(left, bottom),
                     (left, top),
                     (right, top),
                     (right, bottom),
                     (left, bottom)]

            codes = [Path.MOVETO,
                     Path.LINETO,
                     Path.LINETO,
                     Path.LINETO,
                     Path.CLOSEPOLY]

            path = Path(verts, codes)
            ax = fig.add_subplot(111)
            patch = patches.PathPatch(path, fc='white', lw=1)  # , hatch='x')
            ax.add_patch(patch)
        # Plot the channel's cross section from argument, if any
        if x is not None and y is not None:
            plt.plot(x, y)
        ax.axis('equal')
        plt.xlabel(u'Width (m)')
        plt.ylabel(u'Depth (m)')
        plt.title(u'Front view')
        plt.grid(True)
        plt.gca().invert_yaxis()

        # Save fig if there is a path
        if filename == '':
            plt.show()
        else:
            plt.savefig(filename)
            plt.show()
        return True

    def front_view(self, x=None, y=None, filename=''):
        """ Draws the front view of the dam. """

        fig = plt.figure()
        for i in range(len(self.dimensions)):
            left = self.dimensions[i][COORDX]
            right = self.dimensions[i][COORDX] + self.dimensions[i][WIDTH]
            # Because the inverted axis bottom and top should be as follow
            top = self.dimensions[i][COORDY]
            bottom = self.dimensions[i][COORDY] - self.dimensions[i][HEIGHT]

            verts = [(left, bottom),
                     (left, top),
                     (right, top),
                     (right, bottom),
                     (left, bottom)]

            codes = [Path.MOVETO,
                     Path.LINETO,
                     Path.LINETO,
                     Path.LINETO,
                     Path.CLOSEPOLY]

            path = Path(verts, codes)
            ax = fig.add_subplot(111)
            patch = patches.PathPatch(path, fc='white', lw=1)  # , hatch='x')
            ax.add_patch(patch)
        # Plot the channel's cross section from argument, if any
        if x is not None and y is not None:
            plt.plot(x, y)
        ax.axis('equal')
        plt.xlabel(u'Width (m)')
        plt.ylabel(u'Depth (m)')
        plt.title(u'Front view')
        plt.grid(True)
        plt.gca().invert_yaxis()

        # Save fig if there is a path
        if filename == '':
            plt.show()
        else:
            plt.savefig(filename)
            plt.show()
        return True

    def side_view(self, filename=''):
        """ Draws the side view of the dam. """

        fig = plt.figure()
        for i in range(len(self.dimensions)):
            left = 0
            right = self.dimensions[i][LENGTH]
            # Because the inverted axis bottom and top should be as follow
            top = self.dimensions[i][COORDY]
            bottom = self.dimensions[i][COORDY] - self.dimensions[i][HEIGHT]

            verts = [(left, bottom),
                     (left, top),
                     (right, top),
                     (right, bottom),
                     (left, bottom)]

            codes = [Path.MOVETO,
                     Path.LINETO,
                     Path.LINETO,
                     Path.LINETO,
                     Path.CLOSEPOLY]

            path = Path(verts, codes)
            ax = fig.add_subplot(111)
            patch = patches.PathPatch(path, fc='white', lw=1)  # , hatch='x')
            ax.add_patch(patch)
        ax.axis('equal')
        plt.xlabel(u'Length (m)')
        plt.ylabel(u'Depth (m)')
        plt.title(u'Side view')
        plt.grid(True)
        plt.gca().invert_yaxis()

        # Save fig if there is a path
        if filename == '':
            plt.show()
        else:
            plt.savefig(filename)
            plt.show()
        return True

    def top_view(self, filename=''):
        """ Draws the top view of the dam. """

        fig = plt.figure()
        for i in range(len(self.dimensions)):
            left = self.dimensions[i][COORDX]
            right = self.dimensions[i][COORDX] + self.dimensions[i][WIDTH]
            # Because the inverted axis bottom and top should be as follow
            top = 0
            bottom = self.dimensions[i][LENGTH]

            verts = [(left, bottom),
                     (left, top),
                     (right, top),
                     (right, bottom),
                     (left, bottom)]

            codes = [Path.MOVETO,
                     Path.LINETO,
                     Path.LINETO,
                     Path.LINETO,
                     Path.CLOSEPOLY]

            path = Path(verts, codes)
            ax = fig.add_subplot(111)
            patch = patches.PathPatch(path, fc='white', lw=1)  # , hatch='x')
            ax.add_patch(patch)
        ax.axis('equal')
        plt.xlabel(u'Width (m)')
        plt.ylabel(u'Length (m)')
        plt.title(u'Top view')
        plt.grid(True)
        plt.gca().invert_yaxis()

        # Save fig if there is a path
        if filename == '':
            plt.show()
        else:
            plt.savefig(filename)
            plt.show()
        return True

    def report(self):
        """ Report

        Creates a report of the design of the gabion dam
        """
        text = ""
        text += self.channel_report()
        text += self.weir_report()
        text += self.report_dimensions()
        text += self.stability_report()
        return text

# END OF THE CLASSES


if __name__ == "__main__":

    # Choose a layer's height
    heights = [0.15, 0.30, 0.5, 1.0]

    # DATA AND DEFINITIONS

    # Open the file that contain data
    files = ['../examples/data/test_data_01.csv']

    for datafile in files:

        # Rational method
        land_use = 8      # land use (0-13), use show_land_use()
        soil_type = 2     # soil type (0-2), use show_soil_type()
        rainfall = 809.1  # annual precipitation, in mm
        intensity = 111   # rainfall intensity, in mm/h
        wshd_area = 45    # watershed area, in ha
        str_length = 989  # stream length, in m
        str_rise = 23     # stream rise (elevation difference), in m

        # Weir
        weir_len = 2
        weir_c = 1.45
        weir_notch = 'h'
        weir_head = 0  # height of water head, to be calculated

        # OPERATIONS AND CALCULATIONS

        # Datafile
        fdata = DataFile(datafile)
        fdata.read_data()
        # fdata.check_data()

        if not fdata.check_data():
            sys.exit("Input data has wrong format or is corrupt!")

        # Rational method
        rc = RunoffCoefficient(land_use, soil_type, rainfall)
        runoff_coef = rc.get_coefficient()
        # Maximum runoff
        rational = RationalFormula(runoff_coef, intensity, wshd_area)
        ct = ConcentrationTime(wshd_area, str_length, str_rise)
        conc_time = ct.get_concentration_time()
        runoff = rational.get_runoff()

        # Create a gabion dam
        hlayer = heights[3]  # layer's height
        dam = GabionDam(fdata.data, hlayer, runoff, weir_c, weir_len,
                        weir_head, weir_notch)
        if not dam.get_xy():
            sys.exit("Input data has wrong format or is corrupt!")
        dam.get_widths()

        # Create the dimensions of the dam
        dam.design()

        # Perform the stability analysis
        dam.analyze()

        # SHOW RESULTS

        print("\n\n\n")
        print("-" * 40)
        print(" REPORT OF RESULTS")
        print("-" * 40)

        print("Data")
        fdata.show_data()

        print("Rational method")
        print(" Maximum runoff:     {0:.3f} m^3/s".format(runoff))
        print(" Runoff coefficient: {0:.3f}".format(runoff_coef))
        print(" Concentration time: {0:.3f} min.".format(conc_time))

        # dam.show_widths()

        # Create report
        report = dam.report()
        print(report)

        # dam.show_dimensions()
        dam.front_view(dam.data_x, dam.data_y)
        dam.side_view()
        dam.top_view()
