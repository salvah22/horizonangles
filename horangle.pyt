# -*- coding: utf-8 -*-
"""
Pycharm Editor
Author: Salvador Hernandez

This script is intended to address input pre-processing for 'RHESSys' and for other purposes, it calculates,
for a Digital Elevation Model (DEM), the maximum obstruction angular height for each cell in any (polar coordinates)
direction angle, the procedure is reading the DEM as a numpy array (matrix), with pertinent modifications depending on
the direction angle and then calculate the max angular height obstruction for every value in every row, afterwards,
modifying the output as desired, applying a function to the results, scaling the values, in degrees or rads.
"""

import arcpy
import numpy
import os

class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the .pyt file)."""
        self.label = "Horizon Angle Toolbox"
        self.alias = "horangle"

        # List of tool classes associated with this toolbox
        self.tools = [CalculateHorizonAngle]


class CalculateHorizonAngle(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Raster Horizon Angle"
        self.description = "The tool calculates the horizon angle for a DEM in the given direction angle"
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        in_raster = arcpy.Parameter(
            displayName="Input DEM Raster(s)",
            name="in_raster",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input",
            multiValue=True)

        direction_angle = arcpy.Parameter(
            displayName="Direction Angle [0-360]",
            name="direction_angle",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")
        direction_angle.value = "0"

        ws = arcpy.Parameter(
            displayName="Save output to the current ArcGIS workspace",
            name="ws",
            datatype="Boolean",
            parameterType="Optional",
            direction="Input")
        ws.value = True

        # function enabled
        fe = arcpy.Parameter(
            displayName="Apply a function to the output",
            name="fe",
            datatype="Boolean",
            parameterType="Optional",
            direction="Input",
            category="Options",
            enabled=True)
        fe.value = False

        # function operation
        fo = arcpy.Parameter(
            displayName="Function",
            name="fo",
            datatype="GPString",
            parameterType="Optional",
            direction="Input",
            category="Options",
            enabled=False)
        fo.value = "numpy.sin"
        fo.filter.type = "ValueList"
        fo.filter.list = ["numpy.sin", "numpy.cos", "numpy.tan"]

        # scaling enabled
        se = arcpy.Parameter(
            displayName="Apply an operation to the output",
            name="se",
            datatype="Boolean",
            parameterType="Optional",
            direction="Input",
            category="Options",
            enabled=True)
        se.value = False

        # scaling operation
        so = arcpy.Parameter(
            displayName="Operation",
            name="so",
            datatype="GPString",
            parameterType="Optional",
            direction="Input",
            category="Options",
            enabled=False)
        so.value = "multiplication"
        so.filter.type = "ValueList"
        so.filter.list = ["multiplication", "division", "addition", "subtraction"]

        # scaling value
        sv = arcpy.Parameter(
            displayName="Operation value",
            name="sv",
            datatype="GPDouble",
            parameterType="Optional",
            direction="Input",
            category="Options",
            enabled=False)
        sv.value = "1"

        nan_val = arcpy.Parameter(
            displayName="Nodata value (NaN)",
            name="nan_val",
            datatype="GPLong",
            parameterType="Required",
            direction="Input",
            category="Advanced")
        nan_val.value = "-9999"

        # over-writing settings
        ow = arcpy.Parameter(
            displayName="Allow output file overwriting",
            name="ow",
            datatype="Boolean",
            parameterType="Optional",
            direction="Input",
            category="Advanced",
            enabled=True)
        ow.value = False

        degrees = arcpy.Parameter(
            displayName="Result in degrees",
            name="degrees",
            datatype="Boolean",
            parameterType="Optional",
            direction="Input",
            category="Advanced",
            enabled=True)
        degrees.value = True

        # output paths
        out = arcpy.Parameter(
            displayName="Output path(s)",
            name="out",
            datatype="DERasterDataset",
            parameterType="Optional",
            direction="Output",
            category="Advanced",
            multiValue=True)
        out.parameterDependencies = [in_raster.name]
        out.schema.clone = True

        # original workspace
        ows = arcpy.Parameter(
            displayName="original workspace",
            name="ows",
            datatype="GPString",
            parameterType="Derived",
            direction="Output")
        ows.value = arcpy.env.workspace

        parameters = [in_raster, direction_angle, ws, fe, fo, se, so, sv, nan_val, ow, degrees, out, ows]
        return parameters

    def isLicensed(self):  # optional
        # Set whether tool is licensed to execute.
        return True

    def updateParameters(self, parameters):  # optional
        # Modify the values and properties of parameters before internal
        # validation is performed.  This method is called whenever a parameter
        # has been changed.
        """
        this was an experimenet
        if parameters[2].value:  # if current workspace enabled, disable the overwriting param
            parameters[9].enabled = False
            parameters[9].value = False
        else:
            parameters[9].enabled = True
            parameters[9].value = False
        """

        def updateOutputPaths():
            rasterList = parameters[0].valueAsText.split(';')
            outputList = []
            directionAngle = parameters[1].valueAsText

            for raster in rasterList:

                if (raster[0] == raster[-1]) and raster.startswith(("'", '"')):
                    raster = raster[1:-1]
                rasterName, rasterExt = os.path.splitext(os.path.basename(raster))

                if parameters[2].value:  # current workspace
                    arcpy.env.workspace = parameters[12].valueAsText
                    outRasterPath = parameters[12].valueAsText
                else:
                    arcpy.env.workspace = os.path.dirname(raster)
                    outRasterPath = os.path.dirname(raster)

                if outRasterPath[-4:] == ".gdb" or outRasterPath[-5:] == ".gdb" + os.sep:
                    gdb = True
                else:
                    gdb = False

                # rasterName modifications
                if len(rasterExt) == 0 and len(rasterName) >= 8:
                    rasterName = rasterName[:7]

                if len(rasterExt) != 0 and gdb:  # rasterExt != 0
                    outRasterPath = os.path.dirname(outRasterPath)

                # digits
                if len(directionAngle) == 1:
                    outRasterName = rasterName + "_00" + directionAngle + "_"
                elif len(directionAngle) == 2:
                    outRasterName = rasterName + "_0" + directionAngle + "_"
                else:
                    outRasterName = rasterName + "_" + directionAngle + "_"

                # hierarchy for already existing rasters
                n1 = 0
                if not parameters[9].value:
                    while outRasterName + str(n1) + rasterExt in arcpy.ListRasters():
                        n1 += 1
                    while os.path.exists(outRasterPath + os.sep + outRasterName + str(n1) + rasterExt):
                        n1 += 1
                outPath = outRasterPath + os.sep + outRasterName + str(n1) + rasterExt
                outputList.append(outPath)
            arcpy.SetParameterAsText(11, ';'.join(outputList))
            parameters[11].value = ';'.join(outputList)
            return

        if parameters[0].value is None:
            arcpy.SetParameterAsText(11, "")
            parameters[11].value = ""
        else:
            if parameters[0].altered:
                updateOutputPaths()

        if parameters[1].value > 359:
            parameters[1].value = 0
        elif parameters[1].value < 0:
            parameters[1].value = 0

        if parameters[3].value:  # if applying a function is desired
            parameters[4].enabled = True
        else:
            parameters[4].enabled = False
            parameters[4].value = ""

        if parameters[5].value:  # if scalation is desired
            parameters[6].enabled = True
            parameters[7].enabled = True
        else:
            parameters[6].enabled = False
            parameters[6].value = ""
            parameters[7].enabled = False
            parameters[7].value = ""

        arcpy.env.workspace = parameters[12].valueAsText
        return


    def updateMessages(self, parameters):   # optional
        """
        Modify the messages created by internal validation for each tool
        parameter. This method is called after internal validation.
        """
        # parameters[0].setWarningMessage(parameters[0].valueAsText)
        # parameters[11].setWarningMessage(parameters[11].valueAsText)

        # NAN
        if parameters[8].value > 0:
            parameters[8].setWarningMessage('NaN is usually set to negative values')
        # CURRENT WORKSPACE
        if not parameters[2].value:
            parameters[2].setWarningMessage('The output raster will be saved to input raster folder')
        # OVERWRITING
        if parameters[9].value:
            parameters[9].setWarningMessage('Check before running the tool to avoid overwriting your precious files')
        return


    def execute(self, parameters, messages):
        """The source code of the tool."""
        arcpy.env.overwriteOutput = parameters[9].value  # allow overwriting if desired
        arcpy.env.workspace = parameters[12].valueAsText
        if not parameters[9].value:
            messages.addMessage("Overwriting: off")
        rasterList = parameters[0].valueAsText.split(';')
        outPath = parameters[11].valueAsText.split(';')
        messages.addMessage("--------------------------------")

        for rasterNum, raster in enumerate(rasterList):
            # Input/output preparation
            if (raster[0] == raster[-1]) and raster.startswith(("'", '"')):
                raster = raster[1:-1]

            outRasterPath = outPath[rasterNum]
            if (outRasterPath[0] == outRasterPath[-1]) and outRasterPath.startswith(("'", '"')):
                outRasterPath = outRasterPath[1:-1]

            # workspace settings
            arcpy.env.workspace = os.path.dirname(outRasterPath)
            messages.addMessage("--- Workspace for Raster {}/{} ---\n{}".format(rasterNum + 1,
                                                                               len(rasterList), arcpy.env.workspace))
            tempFolder = arcpy.env.workspace
            while tempFolder[-4:] == ".gdb" or tempFolder[-5:] == ".gdb" + os.sep:  # rotation rasters are TIF
                tempFolder = os.path.dirname(tempFolder)
            tempFolder = tempFolder + os.sep + "horangle_temp_folder"  # this is used for rotation rasters
            if not os.path.exists(tempFolder):
                os.makedirs(tempFolder)

            # input raster settings
            inRaster = arcpy.Raster(raster)
            directionAngle = parameters[1].valueAsText
            nanVal = parameters[8].value
            lowerLeft = arcpy.Point(inRaster.extent.XMin, inRaster.extent.YMin)
            cellWidth = inRaster.meanCellWidth
            cellHeight = inRaster.meanCellHeight

            ### Calculating horangle with different approaches depending on the case, angles 0, 90, 180 and 270 usually
            # are used more than intermediate (45) angles and using numpy array manipulations can save time

            if directionAngle == "0" or directionAngle == "90" or directionAngle == "180" or directionAngle == "270":
                rasterArray = arcpy.RasterToNumPyArray(inRaster, nodata_to_value=nanVal)

            """Direction Angle check"""
            if directionAngle == "0":  #east
                cellDistance = cellWidth
            elif directionAngle == "90":  #north
                cellDistance = cellHeight
                rasterArray = numpy.transpose(rasterArray)
            elif directionAngle == "180":  #west
                cellDistance = cellWidth
                rasterArray = numpy.fliplr(rasterArray)
            elif directionAngle == "270":  #south
                cellDistance = cellHeight
                rasterArray = numpy.transpose(rasterArray)
                rasterArray = numpy.fliplr(rasterArray)
            else:
                if 0 < int(directionAngle) < 90 or 180 < int(directionAngle) < 270:
                    cellDistance = cellWidth * numpy.cos(numpy.deg2rad(int(directionAngle) % 90))
                else:
                    cellDistance = cellWidth * numpy.sin(numpy.deg2rad(int(directionAngle) % 90))
                n2 = 0
                while os.path.exists(tempFolder + os.sep + "rotate" + str(n2) + ".tif"):
                    n2 += 1
                tempRasterPath = tempFolder + os.sep + "rotate" + str(n2) + ".tif"
                arcpy.Rotate_management(inRaster, tempRasterPath, directionAngle, "", "CUBIC")
                tempRaster = arcpy.Raster(tempRasterPath)
                rasterArray = arcpy.RasterToNumPyArray(tempRaster, nodata_to_value=nanVal)
                tempLowerLeft = arcpy.Point(tempRaster.extent.XMin, tempRaster.extent.YMin)
                tempCellWidth = tempRaster.meanCellWidth
                try:
                    arcpy.Delete_management(tempRasterPath)
                except:
                    messages.addMessage("Unable to delete the following temporal tif:\n{}".format(tempRasterPath))

            if len(rasterArray.shape) > 2:
                arcpy.AddError("Error: The input raster DEM can only have 1 band")

            outArray = numpy.zeros((rasterArray.shape)) + nanVal  # NAN values array to be filled with horangles

            """Looping over rasterArray rows for angles calculation"""
            n3 = 2
            for row, vals in enumerate(rasterArray):
                for col in range(rasterArray.shape[1] - 1):  # last col doesnt have values to the right
                    """
                    outArray[row, col] = 0  # debugging purposes
                    """
                    # in this line you should correct valsArr to the earth curvature with the according ellipsoid
                    h1 = vals[col]  # h1 is the height in current cell
                    h1_array = numpy.where(vals[col + 1:] != nanVal, vals[col + 1:], 0)  # heights to the right
                    h1_array = h1_array[
                               :numpy.argmax(h1_array) + 1]  # optimization vector crop including h3 (max seen value)
                    seen = set()
                    peaks = []
                    for x, h in enumerate(h1_array):
                        if not h in seen and h > h1:  # we only want the first height appearances
                            peaks.append((h - h1, (x + 1) * cellDistance))
                            seen.add(h)
                    # peaks is a list of tuples, height (dh) and distance (dx) differences for every height>h1
                    if len(seen) == 0:  # h1 is the greatest height looking right
                        outArray[row, col] = 0
                    else:
                        # knowing a height greater than h1 exist, we need to calculate the maximum angular obstruction
                        outArray[row, col] = max([numpy.arctan(dh / dx) for dh, dx in peaks])
                if row == int(rasterArray.shape[0] * n3 / 10):  # {:.2f}
                    messages.addMessage("Raster Array looping row {} ({}%)".format(row, n3 * 10))
                    n3 += 2

            if parameters[3].value:
                if parameters[4].valueAsText == "numpy.sin":
                    outArray = numpy.where(outArray != nanVal, numpy.sin(outArray), nanVal)
                elif parameters[4].valueAsText == "numpy.cos":
                    outArray = numpy.where(outArray != nanVal, numpy.cos(outArray), nanVal)
                elif parameters[4].valueAsText == "numpy.tan":
                    outArray = numpy.where(outArray != nanVal, numpy.tan(outArray), nanVal)

            if parameters[10].value:
                outArray = numpy.where(outArray != nanVal, numpy.rad2deg(outArray), nanVal)

            if parameters[5].value:
                if parameters[6].valueAsText == "multiplication":
                    outArray = numpy.where(outArray != nanVal, outArray * parameters[7].value, nanVal)
                elif parameters[6].valueAsText == "division":
                    outArray = numpy.where(outArray != nanVal, outArray / parameters[7].value, nanVal)
                elif parameters[6].valueAsText == "addition":
                    outArray = numpy.where(outArray != nanVal, outArray + parameters[7].value, nanVal)
                elif parameters[6].valueAsText == "subtraction":
                    outArray = numpy.where(outArray != nanVal, outArray - parameters[7].value, nanVal)

            ###
            if parameters[3].value:
                if parameters[5].value:
                    messages.addMessage("Horizon angle calculation, " + parameters[4].valueAsText + " function & " +
                                        parameters[6].valueAsText + " by " + parameters[7].valueAsText + " OK")
                else:
                    messages.addMessage("Horizon angle calculation and " + parameters[4].valueAsText + " function OK")

            else:
                if parameters[5].value:
                    messages.addMessage("Horizon angle calculation and " + parameters[6].valueAsText + " by " +
                                        parameters[7].valueAsText + " OK")
                else:
                    messages.addMessage("Horizon angle calculation OK")

            if directionAngle == "0":  # east
                pass
            elif directionAngle == "90":  # north
                outArray = numpy.transpose(outArray)
            elif directionAngle == "180":  # west
                outArray = numpy.fliplr(outArray)
            elif directionAngle == "270":  # south
                outArray = numpy.fliplr(outArray)
                outArray = numpy.transpose(outArray)
            else:
                tempRaster2 = arcpy.NumPyArrayToRaster(outArray, tempLowerLeft, tempCellWidth, value_to_nodata=nanVal)
                arcpy.Rotate_management(tempRaster2, outRasterPath, "-" + directionAngle, "", "CUBIC", inRaster.extent)
                # messages.addGPMessages()

            if directionAngle == "0" or directionAngle == "90" or directionAngle == "180" or directionAngle == "270":
                newRaster = arcpy.NumPyArrayToRaster(outArray, lowerLeft, cellWidth, value_to_nodata=nanVal)
                newRaster.save(outRasterPath)


            # arcpy.MakeRasterLayer_management(outRasterPath, os.path.basename(outRasterPath))
            messages.addMessage("Success!\n--- Raster {}/{} created in ---\n{}"
                                " ".format(rasterNum + 1, len(rasterList), outRasterPath))
            messages.addMessage("-----------------------------")

