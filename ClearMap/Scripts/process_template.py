# -*- coding: utf-8 -*-
"""
Template to run the processing pipeline
"""


execfile('./ClearMap/Scripts/parameter_file_template.py')
#
resampleData(**CorrectionResamplingParameterCfos);
resampleData(**CorrectionResamplingParameterAutoFluo);
resampleData(**RegistrationResamplingParameter);
resultDirectory  = alignData(**CorrectionAlignmentParameter);
resultDirectory  = alignData(**RegistrationAlignmentParameter);


#
detectCells(**ImageProcessingParameter);


points, intensities = io.readPoints(ImageProcessingParameter["sink"]);
points, intensities = thresholdPoints(points, intensities, threshold = (20, 900), row = (3,3));
io.writePoints(FilteredCellsFile, (points, intensities));
points = io.readPoints(CorrectionResamplingPointsParameter["pointSource"]);
points = resamplePoints(**CorrectionResamplingPointsParameter);
points = transformPoints(points, transformDirectory = CorrectionAlignmentParameter["resultDirectory"], indices = False, resultDirectory = None);
CorrectionResamplingPointsInverseParameter["pointSource"] = points;
points = resamplePointsInverse(**CorrectionResamplingPointsInverseParameter);
RegistrationResamplingPointParameter["pointSource"] = points;
points = resamplePoints(**RegistrationResamplingPointParameter);
points = transformPoints(points, transformDirectory = RegistrationAlignmentParameter["resultDirectory"], indices = False, resultDirectory = None);
io.writePoints(TransformedCellsFile, points);


# Check Cell; detection
import iDISCO.Visualization.Plot as plt;
pointSource= os.path.join(BaseDirectory, FilteredCellsFile[0]);
data = plt.overlayPoints(cFosFile, pointSource, pointColor = None, **cFosFileRange);
io.writeData(os.path.join(BaseDirectory, 'cells_check.tif'), data);


#################### Heat map generation

points = io.readPoints(TransformedCellsFile)
intensities = io.readPoints(FilteredCellsFile[1])

vox = voxelize(points, AtlasFile, **voxelizeParameter);
io.writeData(os.path.join(BaseDirectory, 'cells_heatmap.tif'), vox.astype('int32'));

voxelizeParameter["weights"] = intensities[:,0].astype(float);
vox = voxelize(points, AtlasFile, **voxelizeParameter);
io.writeData(os.path.join(BaseDirectory, 'cells_heatmap_weighted.tif'), vox.astype('int32'));


ids, counts = countPointsInRegions(points, labeledImage = AnnotationFile, intensities = intensities, intensityRow = 0);
table = numpy.zeros(ids.shape, dtype=[('id','int64'),('counts','f8'),('name', 'a256')])
table["id"] = ids;
table["counts"] = counts;
table["name"] = labelToName(ids);
io.writeTable(os.path.join(BaseDirectory, 'Annotated_counts_intensities.csv'), table);


ids, counts = countPointsInRegions(points, labeledImage = AnnotationFile, intensities = None);
table = numpy.zeros(ids.shape, dtype=[('id','int64'),('counts','f8'),('name', 'a256')])
table["id"] = ids;
table["counts"] = counts;
table["name"] = labelToName(ids);
io.writeTable(os.path.join(BaseDirectory, 'Annotated_counts.csv'), table);



#####################
#####################
#####################
#####################
#####################
#####################
