

minT = 2000;
maxT = 40000


########################
########################

execfile('1whiskers.py')

#resampleData(**CorrectionResamplingParameterCfos);
#resampleData(**CorrectionResamplingParameterAutoFluo);
#resampleData(**RegistrationResamplingParameter);
#resultDirectory  = alignData(**CorrectionAlignmentParameter);
#resultDirectory  = alignData(**RegistrationAlignmentParameter);

detectCells(**ImageProcessingParameter);

points, intensities = io.readPoints((os.path.join(BaseDirectory, 'cells-allpoints.npy'),  os.path.join(BaseDirectory,  'intensities-allpoints.npy')));
points, intensities = thresholdPoints(points, intensities, threshold = (minT, maxT), row = (1,0));
io.writePoints((os.path.join(BaseDirectory, 'cells.npy'), os.path.join(BaseDirectory,  'intensities.npy')), (points, intensities));
points = io.readPoints(CorrectionResamplingPointsParameter["pointSource"]);
points = resamplePoints(**CorrectionResamplingPointsParameter);
points = transformPoints(points, transformDirectory = CorrectionAlignmentParameter["resultDirectory"], indices = False, resultDirectory = None);
CorrectionResamplingPointsInverseParameter["pointSource"] = points;
points = resamplePointsInverse(**CorrectionResamplingPointsInverseParameter);
RegistrationResamplingPointParameter["pointSource"] = points;
points = resamplePoints(**RegistrationResamplingPointParameter);
points = transformPoints(points, transformDirectory = RegistrationAlignmentParameter["resultDirectory"], indices = False, resultDirectory = None);
TransformedCellFile = os.path.join(BaseDirectory, 'cells_transformed_to_Atlas.npy')
io.writePoints(TransformedCellFile, points);

vox = voxelize(points, AtlasFile, **VoxelizationParameter);
io.writeData(os.path.join(BaseDirectory, 'cells_heatmap.tif'), vox.astype('int32'));

VoxelizationParameter["voxelizationWeights"] = intensities[:,0].astype(float);
vox = voxelize(points, AtlasFile, **VoxelizationParameter);
io.writeData(os.path.join(BaseDirectory, 'cells_heatmap_weighted.tif'), vox.astype('int32'));

ids, counts = countPointsInRegions(points, labeledImage = AnnotationFile, intensities = intensities, intensityRow = 0);
table = numpy.zeros(ids.shape, dtype=[('id','int64'),('counts','f8'),('name', 'a256')])
table["id"] = ids;
table["counts"] = counts;
table["name"] = labelToName(ids);
with open(os.path.join(BaseDirectory, 'Annotated_counts_intensities.csv'),'w') as f:
    for sublist in table:
        f.write(', '.join([str(item) for item in sublist]));
        f.write('\n');
    f.close();

ids, counts = countPointsInRegions(points, labeledImage = AnnotationFile, intensities = None);
table = numpy.zeros(ids.shape, dtype=[('id','int64'),('counts','f8'),('name', 'a256')])
table["id"] = ids;
table["counts"] = counts;
table["name"] = labelToName(ids);
with open(os.path.join(BaseDirectory, 'Annotated_counts.csv'),'w') as f:
    for sublist in table:
        f.write(', '.join([str(item) for item in sublist]));
        f.write('\n');
    f.close();


#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################



execfile('2whiskers.py')

#resampleData(**CorrectionResamplingParameterCfos);
#resampleData(**CorrectionResamplingParameterAutoFluo);
#resampleData(**RegistrationResamplingParameter);
#resultDirectory  = alignData(**CorrectionAlignmentParameter);
#resultDirectory  = alignData(**RegistrationAlignmentParameter);

detectCells(**ImageProcessingParameter);

points, intensities = io.readPoints((os.path.join(BaseDirectory, 'cells-allpoints.npy'),  os.path.join(BaseDirectory,  'intensities-allpoints.npy')));
points, intensities = thresholdPoints(points, intensities, threshold = (minT, maxT), row = (1,0));
io.writePoints((os.path.join(BaseDirectory, 'cells.npy'), os.path.join(BaseDirectory,  'intensities.npy')), (points, intensities));
points = io.readPoints(CorrectionResamplingPointsParameter["pointSource"]);
points = resamplePoints(**CorrectionResamplingPointsParameter);
points = transformPoints(points, transformDirectory = CorrectionAlignmentParameter["resultDirectory"], indices = False, resultDirectory = None);
CorrectionResamplingPointsInverseParameter["pointSource"] = points;
points = resamplePointsInverse(**CorrectionResamplingPointsInverseParameter);
RegistrationResamplingPointParameter["pointSource"] = points;
points = resamplePoints(**RegistrationResamplingPointParameter);
points = transformPoints(points, transformDirectory = RegistrationAlignmentParameter["resultDirectory"], indices = False, resultDirectory = None);
TransformedCellFile = os.path.join(BaseDirectory, 'cells_transformed_to_Atlas.npy')
io.writePoints(TransformedCellFile, points);

vox = voxelize(points, AtlasFile, **VoxelizationParameter);
io.writeData(os.path.join(BaseDirectory, 'cells_heatmap.tif'), vox.astype('int32'));

VoxelizationParameter["voxelizationWeights"] = intensities[:,0].astype(float);
vox = voxelize(points, AtlasFile, **VoxelizationParameter);
io.writeData(os.path.join(BaseDirectory, 'cells_heatmap_weighted.tif'), vox.astype('int32'));

ids, counts = countPointsInRegions(points, labeledImage = AnnotationFile, intensities = intensities, intensityRow = 0);
table = numpy.zeros(ids.shape, dtype=[('id','int64'),('counts','f8'),('name', 'a256')])
table["id"] = ids;
table["counts"] = counts;
table["name"] = labelToName(ids);
with open(os.path.join(BaseDirectory, 'Annotated_counts_intensities.csv'),'w') as f:
    for sublist in table:
        f.write(', '.join([str(item) for item in sublist]));
        f.write('\n');
    f.close();

ids, counts = countPointsInRegions(points, labeledImage = AnnotationFile, intensities = None);
table = numpy.zeros(ids.shape, dtype=[('id','int64'),('counts','f8'),('name', 'a256')])
table["id"] = ids;
table["counts"] = counts;
table["name"] = labelToName(ids);
with open(os.path.join(BaseDirectory, 'Annotated_counts.csv'),'w') as f:
    for sublist in table:
        f.write(', '.join([str(item) for item in sublist]));
        f.write('\n');
    f.close();


#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################

execfile('3whiskers.py')

#resampleData(**CorrectionResamplingParameterCfos);
#resampleData(**CorrectionResamplingParameterAutoFluo);
#resampleData(**RegistrationResamplingParameter);
#resultDirectory  = alignData(**CorrectionAlignmentParameter);
#resultDirectory  = alignData(**RegistrationAlignmentParameter);

detectCells(**ImageProcessingParameter);

points, intensities = io.readPoints((os.path.join(BaseDirectory, 'cells-allpoints.npy'),  os.path.join(BaseDirectory,  'intensities-allpoints.npy')));
points, intensities = thresholdPoints(points, intensities, threshold = (minT, maxT), row = (1,0));
io.writePoints((os.path.join(BaseDirectory, 'cells.npy'), os.path.join(BaseDirectory,  'intensities.npy')), (points, intensities));
points = io.readPoints(CorrectionResamplingPointsParameter["pointSource"]);
points = resamplePoints(**CorrectionResamplingPointsParameter);
points = transformPoints(points, transformDirectory = CorrectionAlignmentParameter["resultDirectory"], indices = False, resultDirectory = None);
CorrectionResamplingPointsInverseParameter["pointSource"] = points;
points = resamplePointsInverse(**CorrectionResamplingPointsInverseParameter);
RegistrationResamplingPointParameter["pointSource"] = points;
points = resamplePoints(**RegistrationResamplingPointParameter);
points = transformPoints(points, transformDirectory = RegistrationAlignmentParameter["resultDirectory"], indices = False, resultDirectory = None);
TransformedCellFile = os.path.join(BaseDirectory, 'cells_transformed_to_Atlas.npy')
io.writePoints(TransformedCellFile, points);

vox = voxelize(points, AtlasFile, **VoxelizationParameter);
io.writeData(os.path.join(BaseDirectory, 'cells_heatmap.tif'), vox.astype('int32'));

VoxelizationParameter["voxelizationWeights"] = intensities[:,0].astype(float);
vox = voxelize(points, AtlasFile, **VoxelizationParameter);
io.writeData(os.path.join(BaseDirectory, 'cells_heatmap_weighted.tif'), vox.astype('int32'));

ids, counts = countPointsInRegions(points, labeledImage = AnnotationFile, intensities = intensities, intensityRow = 0);
table = numpy.zeros(ids.shape, dtype=[('id','int64'),('counts','f8'),('name', 'a256')])
table["id"] = ids;
table["counts"] = counts;
table["name"] = labelToName(ids);
with open(os.path.join(BaseDirectory, 'Annotated_counts_intensities.csv'),'w') as f:
    for sublist in table:
        f.write(', '.join([str(item) for item in sublist]));
        f.write('\n');
    f.close();

ids, counts = countPointsInRegions(points, labeledImage = AnnotationFile, intensities = None);
table = numpy.zeros(ids.shape, dtype=[('id','int64'),('counts','f8'),('name', 'a256')])
table["id"] = ids;
table["counts"] = counts;
table["name"] = labelToName(ids);
with open(os.path.join(BaseDirectory, 'Annotated_counts.csv'),'w') as f:
    for sublist in table:
        f.write(', '.join([str(item) for item in sublist]));
        f.write('\n');
    f.close();


#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################

execfile('4whiskers.py')

#resampleData(**CorrectionResamplingParameterCfos);
#resampleData(**CorrectionResamplingParameterAutoFluo);
#resampleData(**RegistrationResamplingParameter);
#resultDirectory  = alignData(**CorrectionAlignmentParameter);
#resultDirectory  = alignData(**RegistrationAlignmentParameter);

detectCells(**ImageProcessingParameter);

points, intensities = io.readPoints((os.path.join(BaseDirectory, 'cells-allpoints.npy'),  os.path.join(BaseDirectory,  'intensities-allpoints.npy')));
points, intensities = thresholdPoints(points, intensities, threshold = (minT, maxT), row = (1,0));
io.writePoints((os.path.join(BaseDirectory, 'cells.npy'), os.path.join(BaseDirectory,  'intensities.npy')), (points, intensities));
points = io.readPoints(CorrectionResamplingPointsParameter["pointSource"]);
points = resamplePoints(**CorrectionResamplingPointsParameter);
points = transformPoints(points, transformDirectory = CorrectionAlignmentParameter["resultDirectory"], indices = False, resultDirectory = None);
CorrectionResamplingPointsInverseParameter["pointSource"] = points;
points = resamplePointsInverse(**CorrectionResamplingPointsInverseParameter);
RegistrationResamplingPointParameter["pointSource"] = points;
points = resamplePoints(**RegistrationResamplingPointParameter);
points = transformPoints(points, transformDirectory = RegistrationAlignmentParameter["resultDirectory"], indices = False, resultDirectory = None);
TransformedCellFile = os.path.join(BaseDirectory, 'cells_transformed_to_Atlas.npy')
io.writePoints(TransformedCellFile, points);

vox = voxelize(points, AtlasFile, **VoxelizationParameter);
io.writeData(os.path.join(BaseDirectory, 'cells_heatmap.tif'), vox.astype('int32'));

VoxelizationParameter["voxelizationWeights"] = intensities[:,0].astype(float);
vox = voxelize(points, AtlasFile, **VoxelizationParameter);
io.writeData(os.path.join(BaseDirectory, 'cells_heatmap_weighted.tif'), vox.astype('int32'));

ids, counts = countPointsInRegions(points, labeledImage = AnnotationFile, intensities = intensities, intensityRow = 0);
table = numpy.zeros(ids.shape, dtype=[('id','int64'),('counts','f8'),('name', 'a256')])
table["id"] = ids;
table["counts"] = counts;
table["name"] = labelToName(ids);
with open(os.path.join(BaseDirectory, 'Annotated_counts_intensities.csv'),'w') as f:
    for sublist in table:
        f.write(', '.join([str(item) for item in sublist]));
        f.write('\n');
    f.close();

ids, counts = countPointsInRegions(points, labeledImage = AnnotationFile, intensities = None);
table = numpy.zeros(ids.shape, dtype=[('id','int64'),('counts','f8'),('name', 'a256')])
table["id"] = ids;
table["counts"] = counts;
table["name"] = labelToName(ids);
with open(os.path.join(BaseDirectory, 'Annotated_counts.csv'),'w') as f:
    for sublist in table:
        f.write(', '.join([str(item) for item in sublist]));
        f.write('\n');
    f.close();


#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################

execfile('5whiskers.py')

#resampleData(**CorrectionResamplingParameterCfos);
#resampleData(**CorrectionResamplingParameterAutoFluo);
#resampleData(**RegistrationResamplingParameter);
#resultDirectory  = alignData(**CorrectionAlignmentParameter);
#resultDirectory  = alignData(**RegistrationAlignmentParameter);

detectCells(**ImageProcessingParameter);

points, intensities = io.readPoints((os.path.join(BaseDirectory, 'cells-allpoints.npy'),  os.path.join(BaseDirectory,  'intensities-allpoints.npy')));
points, intensities = thresholdPoints(points, intensities, threshold = (minT, maxT), row = (1,0));
io.writePoints((os.path.join(BaseDirectory, 'cells.npy'), os.path.join(BaseDirectory,  'intensities.npy')), (points, intensities));
points = io.readPoints(CorrectionResamplingPointsParameter["pointSource"]);
points = resamplePoints(**CorrectionResamplingPointsParameter);
points = transformPoints(points, transformDirectory = CorrectionAlignmentParameter["resultDirectory"], indices = False, resultDirectory = None);
CorrectionResamplingPointsInverseParameter["pointSource"] = points;
points = resamplePointsInverse(**CorrectionResamplingPointsInverseParameter);
RegistrationResamplingPointParameter["pointSource"] = points;
points = resamplePoints(**RegistrationResamplingPointParameter);
points = transformPoints(points, transformDirectory = RegistrationAlignmentParameter["resultDirectory"], indices = False, resultDirectory = None);
TransformedCellFile = os.path.join(BaseDirectory, 'cells_transformed_to_Atlas.npy')
io.writePoints(TransformedCellFile, points);

vox = voxelize(points, AtlasFile, **VoxelizationParameter);
io.writeData(os.path.join(BaseDirectory, 'cells_heatmap.tif'), vox.astype('int32'));

VoxelizationParameter["voxelizationWeights"] = intensities[:,0].astype(float);
vox = voxelize(points, AtlasFile, **VoxelizationParameter);
io.writeData(os.path.join(BaseDirectory, 'cells_heatmap_weighted.tif'), vox.astype('int32'));

ids, counts = countPointsInRegions(points, labeledImage = AnnotationFile, intensities = intensities, intensityRow = 0);
table = numpy.zeros(ids.shape, dtype=[('id','int64'),('counts','f8'),('name', 'a256')])
table["id"] = ids;
table["counts"] = counts;
table["name"] = labelToName(ids);
with open(os.path.join(BaseDirectory, 'Annotated_counts_intensities.csv'),'w') as f:
    for sublist in table:
        f.write(', '.join([str(item) for item in sublist]));
        f.write('\n');
    f.close();

ids, counts = countPointsInRegions(points, labeledImage = AnnotationFile, intensities = None);
table = numpy.zeros(ids.shape, dtype=[('id','int64'),('counts','f8'),('name', 'a256')])
table["id"] = ids;
table["counts"] = counts;
table["name"] = labelToName(ids);
with open(os.path.join(BaseDirectory, 'Annotated_counts.csv'),'w') as f:
    for sublist in table:
        f.write(', '.join([str(item) for item in sublist]));
        f.write('\n');
    f.close();


#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################

execfile('6whiskers.py')

#resampleData(**CorrectionResamplingParameterCfos);
#resampleData(**CorrectionResamplingParameterAutoFluo);
#resampleData(**RegistrationResamplingParameter);
#resultDirectory  = alignData(**CorrectionAlignmentParameter);
#resultDirectory  = alignData(**RegistrationAlignmentParameter);

#detectCells(**ImageProcessingParameter);

points, intensities = io.readPoints((os.path.join(BaseDirectory, 'cells-allpoints.npy'),  os.path.join(BaseDirectory,  'intensities-allpoints.npy')));
points, intensities = thresholdPoints(points, intensities, threshold = (minT, maxT), row = (1,0));
io.writePoints((os.path.join(BaseDirectory, 'cells.npy'), os.path.join(BaseDirectory,  'intensities.npy')), (points, intensities));
points = io.readPoints(CorrectionResamplingPointsParameter["pointSource"]);
points = resamplePoints(**CorrectionResamplingPointsParameter);
points = transformPoints(points, transformDirectory = CorrectionAlignmentParameter["resultDirectory"], indices = False, resultDirectory = None);
CorrectionResamplingPointsInverseParameter["pointSource"] = points;
points = resamplePointsInverse(**CorrectionResamplingPointsInverseParameter);
RegistrationResamplingPointParameter["pointSource"] = points;
points = resamplePoints(**RegistrationResamplingPointParameter);
points = transformPoints(points, transformDirectory = RegistrationAlignmentParameter["resultDirectory"], indices = False, resultDirectory = None);
TransformedCellFile = os.path.join(BaseDirectory, 'cells_transformed_to_Atlas.npy')
io.writePoints(TransformedCellFile, points);

vox = voxelize(points, AtlasFile, **VoxelizationParameter);
io.writeData(os.path.join(BaseDirectory, 'cells_heatmap.tif'), vox.astype('int32'));

VoxelizationParameter["voxelizationWeights"] = intensities[:,0].astype(float);
vox = voxelize(points, AtlasFile, **VoxelizationParameter);
io.writeData(os.path.join(BaseDirectory, 'cells_heatmap_weighted.tif'), vox.astype('int32'));

ids, counts = countPointsInRegions(points, labeledImage = AnnotationFile, intensities = intensities, intensityRow = 0);
table = numpy.zeros(ids.shape, dtype=[('id','int64'),('counts','f8'),('name', 'a256')])
table["id"] = ids;
table["counts"] = counts;
table["name"] = labelToName(ids);
with open(os.path.join(BaseDirectory, 'Annotated_counts_intensities.csv'),'w') as f:
    for sublist in table:
        f.write(', '.join([str(item) for item in sublist]));
        f.write('\n');
    f.close();

ids, counts = countPointsInRegions(points, labeledImage = AnnotationFile, intensities = None);
table = numpy.zeros(ids.shape, dtype=[('id','int64'),('counts','f8'),('name', 'a256')])
table["id"] = ids;
table["counts"] = counts;
table["name"] = labelToName(ids);
with open(os.path.join(BaseDirectory, 'Annotated_counts.csv'),'w') as f:
    for sublist in table:
        f.write(', '.join([str(item) for item in sublist]));
        f.write('\n');
    f.close();


#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################

execfile('7whiskers.py')

#resampleData(**CorrectionResamplingParameterCfos);
#resampleData(**CorrectionResamplingParameterAutoFluo);
#resampleData(**RegistrationResamplingParameter);
#resultDirectory  = alignData(**CorrectionAlignmentParameter);
#resultDirectory  = alignData(**RegistrationAlignmentParameter);

#detectCells(**ImageProcessingParameter);

points, intensities = io.readPoints((os.path.join(BaseDirectory, 'cells-allpoints.npy'),  os.path.join(BaseDirectory,  'intensities-allpoints.npy')));
points, intensities = thresholdPoints(points, intensities, threshold = (minT, maxT), row = (1,0));
io.writePoints((os.path.join(BaseDirectory, 'cells.npy'), os.path.join(BaseDirectory,  'intensities.npy')), (points, intensities));
points = io.readPoints(CorrectionResamplingPointsParameter["pointSource"]);
points = resamplePoints(**CorrectionResamplingPointsParameter);
points = transformPoints(points, transformDirectory = CorrectionAlignmentParameter["resultDirectory"], indices = False, resultDirectory = None);
CorrectionResamplingPointsInverseParameter["pointSource"] = points;
points = resamplePointsInverse(**CorrectionResamplingPointsInverseParameter);
RegistrationResamplingPointParameter["pointSource"] = points;
points = resamplePoints(**RegistrationResamplingPointParameter);
points = transformPoints(points, transformDirectory = RegistrationAlignmentParameter["resultDirectory"], indices = False, resultDirectory = None);
TransformedCellFile = os.path.join(BaseDirectory, 'cells_transformed_to_Atlas.npy')
io.writePoints(TransformedCellFile, points);

vox = voxelize(points, AtlasFile, **VoxelizationParameter);
io.writeData(os.path.join(BaseDirectory, 'cells_heatmap.tif'), vox.astype('int32'));

VoxelizationParameter["voxelizationWeights"] = intensities[:,0].astype(float);
vox = voxelize(points, AtlasFile, **VoxelizationParameter);
io.writeData(os.path.join(BaseDirectory, 'cells_heatmap_weighted.tif'), vox.astype('int32'));

ids, counts = countPointsInRegions(points, labeledImage = AnnotationFile, intensities = intensities, intensityRow = 0);
table = numpy.zeros(ids.shape, dtype=[('id','int64'),('counts','f8'),('name', 'a256')])
table["id"] = ids;
table["counts"] = counts;
table["name"] = labelToName(ids);
with open(os.path.join(BaseDirectory, 'Annotated_counts_intensities.csv'),'w') as f:
    for sublist in table:
        f.write(', '.join([str(item) for item in sublist]));
        f.write('\n');
    f.close();

ids, counts = countPointsInRegions(points, labeledImage = AnnotationFile, intensities = None);
table = numpy.zeros(ids.shape, dtype=[('id','int64'),('counts','f8'),('name', 'a256')])
table["id"] = ids;
table["counts"] = counts;
table["name"] = labelToName(ids);
with open(os.path.join(BaseDirectory, 'Annotated_counts.csv'),'w') as f:
    for sublist in table:
        f.write(', '.join([str(item) for item in sublist]));
        f.write('\n');
    f.close();


#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################

execfile('9whiskers.py')

#resampleData(**CorrectionResamplingParameterCfos);
#resampleData(**CorrectionResamplingParameterAutoFluo);
#resampleData(**RegistrationResamplingParameter);
#resultDirectory  = alignData(**CorrectionAlignmentParameter);
#resultDirectory  = alignData(**RegistrationAlignmentParameter);

#detectCells(**ImageProcessingParameter);


points, intensities = io.readPoints((os.path.join(BaseDirectory, 'cells-allpoints.npy'),  os.path.join(BaseDirectory,  'intensities-allpoints.npy')));
points, intensities = thresholdPoints(points, intensities, threshold = (minT, maxT), row = (1,0));
io.writePoints((os.path.join(BaseDirectory, 'cells.npy'), os.path.join(BaseDirectory,  'intensities.npy')), (points, intensities));
points = io.readPoints(CorrectionResamplingPointsParameter["pointSource"]);
points = resamplePoints(**CorrectionResamplingPointsParameter);
points = transformPoints(points, transformDirectory = CorrectionAlignmentParameter["resultDirectory"], indices = False, resultDirectory = None);
CorrectionResamplingPointsInverseParameter["pointSource"] = points;
points = resamplePointsInverse(**CorrectionResamplingPointsInverseParameter);
RegistrationResamplingPointParameter["pointSource"] = points;
points = resamplePoints(**RegistrationResamplingPointParameter);
points = transformPoints(points, transformDirectory = RegistrationAlignmentParameter["resultDirectory"], indices = False, resultDirectory = None);
TransformedCellFile = os.path.join(BaseDirectory, 'cells_transformed_to_Atlas.npy')
io.writePoints(TransformedCellFile, points);

vox = voxelize(points, AtlasFile, **VoxelizationParameter);
io.writeData(os.path.join(BaseDirectory, 'cells_heatmap.tif'), vox.astype('int32'));

VoxelizationParameter["voxelizationWeights"] = intensities[:,0].astype(float);
vox = voxelize(points, AtlasFile, **VoxelizationParameter);
io.writeData(os.path.join(BaseDirectory, 'cells_heatmap_weighted.tif'), vox.astype('int32'));

ids, counts, countsi = countPointsInRegions(points, labeledImage = AnnotationFile, intensities = intensities, intensityRow = 0, returnCounts = True);
table = numpy.zeros(ids.shape, dtype=[('id','int64'),('counts','f8'),('name', 'a256')])
table["id"] = ids;
table["counts"] = countsi;
table["name"] = labelToName(ids);
with open(os.path.join(BaseDirectory, 'Annotated_counts_intensities.csv'),'w') as f:
    for sublist in table:
        f.write(', '.join([str(item) for item in sublist]));
        f.write('\n');
    f.close();

#ids, counts = countPointsInRegions(points, labeledImage = AnnotationFile, intensities = None);
table = numpy.zeros(ids.shape, dtype=[('id','int64'),('counts','f8'),('name', 'a256')])
table["id"] = ids;
table["counts"] = counts;
table["name"] = labelToName(ids);
with open(os.path.join(BaseDirectory, 'Annotated_counts.csv'),'w') as f:
    for sublist in table:
        f.write(', '.join([str(item) for item in sublist]));
        f.write('\n');
    f.close();


#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################

execfile('10whiskers.py')

#resampleData(**CorrectionResamplingParameterCfos);
#resampleData(**CorrectionResamplingParameterAutoFluo);
#resampleData(**RegistrationResamplingParameter);
#resultDirectory  = alignData(**CorrectionAlignmentParameter);
#resultDirectory  = alignData(**RegistrationAlignmentParameter);

detectCells(**ImageProcessingParameter);

points, intensities = io.readPoints((os.path.join(BaseDirectory, 'cells-allpoints.npy'),  os.path.join(BaseDirectory,  'intensities-allpoints.npy')));
points, intensities = thresholdPoints(points, intensities, threshold = (minT, maxT), row = (1,0));
io.writePoints((os.path.join(BaseDirectory, 'cells.npy'), os.path.join(BaseDirectory,  'intensities.npy')), (points, intensities));
points = io.readPoints(CorrectionResamplingPointsParameter["pointSource"]);
points = resamplePoints(**CorrectionResamplingPointsParameter);
points = transformPoints(points, transformDirectory = CorrectionAlignmentParameter["resultDirectory"], indices = False, resultDirectory = None);
CorrectionResamplingPointsInverseParameter["pointSource"] = points;
points = resamplePointsInverse(**CorrectionResamplingPointsInverseParameter);
RegistrationResamplingPointParameter["pointSource"] = points;
points = resamplePoints(**RegistrationResamplingPointParameter);
points = transformPoints(points, transformDirectory = RegistrationAlignmentParameter["resultDirectory"], indices = False, resultDirectory = None);
TransformedCellFile = os.path.join(BaseDirectory, 'cells_transformed_to_Atlas.npy')
io.writePoints(TransformedCellFile, points);

vox = voxelize(points, AtlasFile, **VoxelizationParameter);
io.writeData(os.path.join(BaseDirectory, 'cells_heatmap.tif'), vox.astype('int32'));

VoxelizationParameter["voxelizationWeights"] = intensities[:,0].astype(float);
vox = voxelize(points, AtlasFile, **VoxelizationParameter);
io.writeData(os.path.join(BaseDirectory, 'cells_heatmap_weighted.tif'), vox.astype('int32'));

ids, counts = countPointsInRegions(points, labeledImage = AnnotationFile, intensities = intensities, intensityRow = 0);
table = numpy.zeros(ids.shape, dtype=[('id','int64'),('counts','f8'),('name', 'a256')])
table["id"] = ids;
table["counts"] = counts;
table["name"] = labelToName(ids);
with open(os.path.join(BaseDirectory, 'Annotated_counts_intensities.csv'),'w') as f:
    for sublist in table:
        f.write(', '.join([str(item) for item in sublist]));
        f.write('\n');
    f.close();

ids, counts = countPointsInRegions(points, labeledImage = AnnotationFile, intensities = None);
table = numpy.zeros(ids.shape, dtype=[('id','int64'),('counts','f8'),('name', 'a256')])
table["id"] = ids;
table["counts"] = counts;
table["name"] = labelToName(ids);
with open(os.path.join(BaseDirectory, 'Annotated_counts.csv'),'w') as f:
    for sublist in table:
        f.write(', '.join([str(item) for item in sublist]));
        f.write('\n');
    f.close();


#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################