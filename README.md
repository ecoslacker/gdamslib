# gdamslib
Design of small gabion dams for soil conservation

**WARNING!** This module is under *testing*, please consider that there may be calculation errors.
Do not use under production or academic environments until the code is marked as *stable*.

# Description

**gdamslib** is a Pythom module to design small gabion dams for soil conservation purposes.
It can be used in programs that run in command line interface and can be linked to programs
with graphic interface.

The module has modules to perform the following calculations:
* Runoff estimation using the Rational Formula.
* Runoff estimation using Chezy-Manning equation.
* Design of a broad crested rectangular weir.
* Dam dimensioning.
* Stability analysis using the method of Oropeza-Mota & Lopez-Martinez (2009).
* Show graphics of the dam geometry.

# Requirements

* Python 2.7

# Usage

The following code presents an exampleof how to use the module to design one or mode gabion dams
usin the data from a list of file names.

```python
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
```

# License


Design of small gabion dams for soil conservation

Copyright (c) 2014-2018 Eduardo Jiménez &lt;ecoslacker@irriapps.com&gt;

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
SOFTWARE`

