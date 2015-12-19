Search.setIndex({envversion:46,filenames:["api/ClearMap","api/ClearMap.Alignment","api/ClearMap.Analysis","api/ClearMap.Analysis.Tools","api/ClearMap.IO","api/ClearMap.ImageProcessing","api/ClearMap.ImageProcessing.Filter","api/ClearMap.Parameter","api/ClearMap.Settings","api/ClearMap.Utils","api/ClearMap.Visualization","imageanalysis","index","introduction"],objects:{"":{ClearMap:[0,0,0,"-"]},"ClearMap.Alignment":{Elastix:[1,0,0,"-"],Resampling:[1,0,0,"-"]},"ClearMap.Alignment.Elastix":{ElastixBinary:[1,3,1,""],ElastixLib:[1,3,1,""],Initialized:[1,3,1,""],TransformixBinary:[1,3,1,""],alignData:[1,1,1,""],checkElastixInitialized:[1,1,1,""],deformationDistance:[1,1,1,""],deformationField:[1,1,1,""],getResultDataFile:[1,1,1,""],getTransformFileSizeAndSpacing:[1,1,1,""],getTransformParameterFile:[1,1,1,""],initializeElastix:[1,1,1,""],parseElastixOutputPoints:[1,1,1,""],printSettings:[1,1,1,""],rescaleSizeAndSpacing:[1,1,1,""],setElastixLibraryPath:[1,1,1,""],setPathTransformParameterFiles:[1,1,1,""],setTransformFileSizeAndSpacing:[1,1,1,""],test:[1,1,1,""],transformData:[1,1,1,""],transformPoints:[1,1,1,""],writePoints:[1,1,1,""]},"ClearMap.Alignment.Resampling":{fixInterpolation:[1,1,1,""],fixOrientation:[1,1,1,""],inverseOrientation:[1,1,1,""],orientDataSize:[1,1,1,""],orientDataSizeInverse:[1,1,1,""],orientResolution:[1,1,1,""],orientResolutionInverse:[1,1,1,""],orientationToPermuation:[1,1,1,""],resampleData:[1,1,1,""],resampleDataInverse:[1,1,1,""],resampleDataSize:[1,1,1,""],resamplePoints:[1,1,1,""],resamplePointsInverse:[1,1,1,""],resampleXY:[1,1,1,""],sagittalToCoronalData:[1,1,1,""]},"ClearMap.Analysis":{Label:[2,0,0,"-"],Statistics:[2,0,0,"-"],Tools:[3,0,0,"-"],Voxelization:[2,0,0,"-"]},"ClearMap.Analysis.Label":{DefaultAnnotationFile:[2,3,1,""],DefaultLabeledImageFile:[2,3,1,""],Label:[2,3,1,""],LabelInfo:[2,5,1,""],LabelRecord:[2,5,1,""],countPointsInRegions:[2,1,1,""],initialize:[2,1,1,""],labelAtCollapse:[2,1,1,""],labelAtLevel:[2,1,1,""],labelPoints:[2,1,1,""],labelToAcronym:[2,1,1,""],labelToColor:[2,1,1,""],labelToName:[2,1,1,""],makeColorAnnotations:[2,1,1,""],makeColorPalette:[2,1,1,""],test:[2,1,1,""],writeLUT:[2,1,1,""],writePAL:[2,1,1,""]},"ClearMap.Analysis.Label.LabelInfo":{acronym:[2,2,1,""],acronyms:[2,4,1,""],collapse:[2,4,1,""],collapseMap:[2,4,1,""],color:[2,2,1,""],colors:[2,4,1,""],ids:[2,4,1,""],initialize:[2,2,1,""],level:[2,2,1,""],levels:[2,4,1,""],name:[2,2,1,""],names:[2,4,1,""],parent:[2,2,1,""],parents:[2,4,1,""],toLabelAtCollapse:[2,2,1,""],toLabelAtCollapseMap:[2,2,1,""],toLabelAtLevel:[2,2,1,""]},"ClearMap.Analysis.Label.LabelRecord":{"__getnewargs__":[2,2,1,""],"__getstate__":[2,2,1,""],"__repr__":[2,2,1,""],acronym:[2,4,1,""],collapse:[2,4,1,""],color:[2,4,1,""],id:[2,4,1,""],name:[2,4,1,""],parent:[2,4,1,""]},"ClearMap.Analysis.Statistics":{"var":[2,1,1,""],colorPValues:[2,1,1,""],countPointsGroupInRegions:[2,1,1,""],cutoffPValues:[2,1,1,""],mean:[2,1,1,""],readDataGroup:[2,1,1,""],readPointsGroup:[2,1,1,""],std:[2,1,1,""],tTestPointsInRegions:[2,1,1,""],tTestVoxelization:[2,1,1,""],test:[2,1,1,""],testCompletedCumulatives:[2,1,1,""],testCompletedCumulativesInSpheres:[2,1,1,""],testCompletedInvertedCumulatives:[2,1,1,""],thresholdPoints:[2,1,1,""],weightsFromPrecentiles:[2,1,1,""]},"ClearMap.Analysis.Tools":{Extrapolate:[3,0,0,"-"],MultipleComparisonCorrection:[3,0,0,"-"],StatisticalTests:[3,0,0,"-"]},"ClearMap.Analysis.Tools.Extrapolate":{extrap1d:[3,1,1,""],extrap1dFromInterp1d:[3,1,1,""]},"ClearMap.Analysis.Tools.MultipleComparisonCorrection":{correctPValues:[3,1,1,""],estimateQValues:[3,1,1,""]},"ClearMap.Analysis.Tools.StatisticalTests":{testCramerVonMises2Sample:[3,1,1,""]},"ClearMap.Analysis.Voxelization":{test:[2,1,1,""],voxelize:[2,1,1,""],voxelizePixel:[2,1,1,""]},"ClearMap.IO":{CSV:[4,0,0,"-"],FileList:[4,0,0,"-"],IO:[4,0,0,"-"],Imaris:[4,0,0,"-"],NPY:[4,0,0,"-"],NRRD:[4,0,0,"-"],RAW:[4,0,0,"-"],TIF:[4,0,0,"-"],VTK:[4,0,0,"-"]},"ClearMap.IO.CSV":{readPoints:[4,1,1,""],test:[4,1,1,""],writePoints:[4,1,1,""]},"ClearMap.IO.FileList":{copyData:[4,1,1,""],dataSize:[4,1,1,""],dataZSize:[4,1,1,""],fileExperssionToFileName:[4,1,1,""],readData:[4,1,1,""],readDataFiles:[4,1,1,""],readFileList:[4,1,1,""],splitFileExpression:[4,1,1,""],test:[4,1,1,""],writeData:[4,1,1,""]},"ClearMap.IO.IO":{convertData:[4,1,1,""],copyData:[4,1,1,""],copyFile:[4,1,1,""],createDirectory:[4,1,1,""],dataFileExtensionToType:[4,3,1,""],dataFileExtensions:[4,3,1,""],dataFileNameToModule:[4,1,1,""],dataFileNameToType:[4,1,1,""],dataFileTypes:[4,3,1,""],dataSize:[4,1,1,""],dataSizeFromDataRange:[4,1,1,""],dataToRange:[4,1,1,""],dataZSize:[4,1,1,""],fileExtension:[4,1,1,""],isDataFile:[4,1,1,""],isFile:[4,1,1,""],isFileExpression:[4,1,1,""],isPointFile:[4,1,1,""],pointFileExtensionToType:[4,3,1,""],pointFileExtensions:[4,3,1,""],pointFileNameToModule:[4,1,1,""],pointFileNameToType:[4,1,1,""],pointFileTypes:[4,3,1,""],pointShiftFromRange:[4,1,1,""],pointsToCoordinates:[4,1,1,""],pointsToCoordinatesAndProperties:[4,1,1,""],pointsToCoordinatesAndPropertiesFileNames:[4,1,1,""],pointsToProperties:[4,1,1,""],pointsToRange:[4,1,1,""],readData:[4,1,1,""],readPoints:[4,1,1,""],toDataRange:[4,1,1,""],toDataSize:[4,1,1,""],toMultiChannelData:[4,1,1,""],writeData:[4,1,1,""],writePoints:[4,1,1,""],writeTable:[4,1,1,""]},"ClearMap.IO.Imaris":{closeFile:[4,1,1,""],copyData:[4,1,1,""],dataSize:[4,1,1,""],dataZSize:[4,1,1,""],getDataExtent:[4,1,1,""],getDataSize:[4,1,1,""],getScaleAndOffset:[4,1,1,""],openFile:[4,1,1,""],readData:[4,1,1,""],readDataSet:[4,1,1,""],readPoints:[4,1,1,""],test:[4,1,1,""],transformPointsToImaris:[4,1,1,""],writeData:[4,1,1,""],writePoints:[4,1,1,""]},"ClearMap.IO.NPY":{readPoints:[4,1,1,""],test:[4,1,1,""],writePoints:[4,1,1,""]},"ClearMap.IO.NRRD":{NrrdError:[4,6,1,""],copyData:[4,1,1,""],dataSize:[4,1,1,""],dataZSize:[4,1,1,""],parse_nrrdvector:[4,1,1,""],parse_optional_nrrdvector:[4,1,1,""],readData:[4,1,1,""],readHeader:[4,1,1,""],test:[4,1,1,""],writeData:[4,1,1,""]},"ClearMap.IO.RAW":{copyData:[4,1,1,""],dataSize:[4,1,1,""],dataZSize:[4,1,1,""],readData:[4,1,1,""],test:[4,1,1,""],writeData:[4,1,1,""],writeHeader:[4,1,1,""],writeRawData:[4,1,1,""]},"ClearMap.IO.TIF":{copyData:[4,1,1,""],dataSize:[4,1,1,""],dataZSize:[4,1,1,""],readData:[4,1,1,""],test:[4,1,1,""],writeData:[4,1,1,""]},"ClearMap.IO.VTK":{readPoints:[4,1,1,""],writePoints:[4,1,1,""]},"ClearMap.ImageProcessing":{BackgroundRemoval:[5,0,0,"-"],CellDetection:[5,0,0,"-"],CellSizeDetection:[5,0,0,"-"],Filter:[6,0,0,"-"],GreyReconstruction:[5,0,0,"-"],IlastikClassification:[5,0,0,"-"],IlluminationCorrection:[5,0,0,"-"],ImageStatistics:[5,0,0,"-"],MaximaDetection:[5,0,0,"-"],SpotDetection:[5,0,0,"-"],StackProcessing:[5,0,0,"-"]},"ClearMap.ImageProcessing.BackgroundRemoval":{removeBackground:[5,1,1,""]},"ClearMap.ImageProcessing.CellDetection":{detectCells:[5,1,1,""]},"ClearMap.ImageProcessing.CellSizeDetection":{detectCellShape:[5,1,1,""],findCellIntensity:[5,1,1,""],findCellSize:[5,1,1,""]},"ClearMap.ImageProcessing.Filter":{Convolution:[6,0,0,"-"],DoGFilter:[6,0,0,"-"],FilterKernel:[6,0,0,"-"],LinearFilter:[6,0,0,"-"],StructureElement:[6,0,0,"-"]},"ClearMap.ImageProcessing.Filter.Convolution":{convolve:[6,1,1,""]},"ClearMap.ImageProcessing.Filter.DoGFilter":{filterDoG:[6,1,1,""]},"ClearMap.ImageProcessing.Filter.FilterKernel":{filterKernel2D:[6,1,1,""],filterKernel3D:[6,1,1,""],filterKernel:[6,1,1,""],test:[6,1,1,""]},"ClearMap.ImageProcessing.Filter.LinearFilter":{filterLinear:[6,1,1,""]},"ClearMap.ImageProcessing.Filter.StructureElement":{structureElement2D:[6,1,1,""],structureElement3D:[6,1,1,""],structureElement:[6,1,1,""],structureElementOffsets:[6,1,1,""]},"ClearMap.ImageProcessing.GreyReconstruction":{greyReconstruction:[5,1,1,""],reconstruct:[5,1,1,""]},"ClearMap.ImageProcessing.IlastikClassification":{Initialized:[5,3,1,""],classifyCells:[5,1,1,""],classifyPixel:[5,1,1,""],initializeIlastik:[5,1,1,""],rescaleFactorIlastik:[5,1,1,""],rescaleToIlastik:[5,1,1,""]},"ClearMap.ImageProcessing.IlluminationCorrection":{DefaultFlatFieldLineFile:[5,3,1,""],correctIllumination:[5,1,1,""],flatfieldFromLine:[5,1,1,""],flatfieldLineFromRegression:[5,1,1,""]},"ClearMap.ImageProcessing.ImageStatistics":{calculateStatistics:[5,1,1,""],calculateStatisticsOnStack:[5,1,1,""],joinStatistics:[5,1,1,""]},"ClearMap.ImageProcessing.MaximaDetection":{extendedMax:[5,1,1,""],findCenterOfMaxima:[5,1,1,""],findExtendedMaxima:[5,1,1,""],findIntensity:[5,1,1,""],findPixelCoordinates:[5,1,1,""],hMaxTransform:[5,1,1,""],localMax:[5,1,1,""]},"ClearMap.ImageProcessing.SpotDetection":{detectSpots:[5,1,1,""],test:[5,1,1,""]},"ClearMap.ImageProcessing.StackProcessing":{calculateChunkSize:[5,1,1,""],calculateSubStacks:[5,1,1,""],joinPoints:[5,1,1,""],noProcessing:[5,1,1,""],parallelProcessStack:[5,1,1,""],printSubStackInfo:[5,1,1,""],sequentiallyProcessStack:[5,1,1,""],writeSubStack:[5,1,1,""]},"ClearMap.Parameter":{AlignmentParameter:[7,3,1,""],IlastikParameter:[7,3,1,""],ResamplingParameter:[7,3,1,""],VoxelizationParameter:[7,3,1,""],detectCellParameter:[7,3,1,""],processStackParameter:[7,3,1,""]},"ClearMap.Settings":{ClearMapPath:[8,3,1,""],ElastixPath:[8,3,1,""],IlastikPath:[8,3,1,""],clearMapPath:[8,1,1,""],setup:[8,1,1,""]},"ClearMap.Utils":{ParameterTools:[9,0,0,"-"],ProcessWriter:[9,0,0,"-"],Timer:[9,0,0,"-"]},"ClearMap.Utils.ParameterTools":{getParameter:[9,1,1,""],joinParameter:[9,1,1,""],writeParameter:[9,1,1,""]},"ClearMap.Utils.ProcessWriter":{ProcessWriter:[9,5,1,""]},"ClearMap.Utils.ProcessWriter.ProcessWriter":{process:[9,4,1,""],write:[9,2,1,""],writeString:[9,2,1,""]},"ClearMap.Utils.Timer":{Timer:[9,5,1,""]},"ClearMap.Utils.Timer.Timer":{elapsedTime:[9,2,1,""],formatElapsedTime:[9,2,1,""],printElapsedTime:[9,2,1,""],reset:[9,2,1,""],start:[9,2,1,""],time:[9,4,1,""]},"ClearMap.Visualization":{Plot:[10,0,0,"-"]},"ClearMap.Visualization.Plot":{overlayLabel:[10,1,1,""],overlayPoints:[10,1,1,""],plotOverlayLabel:[10,1,1,""],plotOverlayPoints:[10,1,1,""],plotTiling:[10,1,1,""],test:[10,1,1,""]},ClearMap:{Alignment:[1,0,0,"-"],Analysis:[2,0,0,"-"],IO:[4,0,0,"-"],ImageProcessing:[5,0,0,"-"],Parameter:[7,0,0,"-"],Settings:[8,0,0,"-"],Utils:[9,0,0,"-"],Visualization:[10,0,0,"-"]}},objnames:{"0":["py","module","Python module"],"1":["py","function","Python function"],"2":["py","method","Python method"],"3":["py","data","Python data"],"4":["py","attribute","Python attribute"],"5":["py","class","Python class"],"6":["py","exception","Python exception"]},objtypes:{"0":"py:module","1":"py:function","2":"py:method","3":"py:data","4":"py:attribute","5":"py:class","6":"py:exception"},terms:{"0x7fee9c25dd70":11,"150524_0_8x":7,"20hf_ultraii_c00_xyz":1,"20hfautofluor_18":7,"25um":2,"27_0_8x":1,"2nd":5,"__getnewargs__":2,"__getstate__":2,"__repr__":2,"_intens":4,"case":1,"class":[2,4,5,9],"const":1,"default":[1,2,3,5,7,9],"final":[2,5,11],"float":[1,3,4,5,6,9,10,11],"function":[1,3,4,5,7,9,10,11,12],"import":[1,4,5,11,12],"int":[1,3,4,5,6,9,11],"new":[0,1,4,5,6,9,10,12],"null":3,"public":[0,12],"return":[1,2,3,4,5,6,8,9,10],"static":5,"true":[1,2,3,4,5,7,9,11],"try":5,"var":2,"while":[1,5],aba:13,abolsut:1,aboslut:1,about:[5,11],absolu:8,absolut:[4,5,8],accept:4,acces:4,access:12,accord:1,accordingli:5,account:1,accuar:11,achiev:5,acronym:2,activ:[6,13],actual:[4,11],adapt:5,add:[1,4],addit:[4,5,13],adjust:2,adult:13,advis:5,affin:[1,7],affineparameterfil:[1,7],after:[4,5],alan:[0,12],algorithm:[5,7],alia:2,align:0,aligndata:1,alignmentdirectori:7,alignmentparamet:7,all:[1,2,3,4,5,7,10,11],allan:13,allen:2,allid:2,allow:[2,5,13],along:5,alpha:10,alreadi:4,also:[4,5,11,13],among:3,analysi:0,analysis:[2,13],analyz:13,andersondarl:2,andresondarl:2,ani:[4,5,6,13],annoat:[2,13],annot:[2,4,13],annotatd:2,annotation_25_right:2,annotationfil:2,anot:[2,13],antibodi:13,antigen:13,api:13,appli:[1,6],applic:[5,6],ara2_annotation_info:2,ara2_annotation_info_collaps:2,arang:[4,13],area:2,arg:[1,2,4,5,9],argument:[4,9],around:5,arrai:[1,2,3,4,5,6,10,13],arrang:1,arrrai:4,assai:13,assembla:12,assing:2,asstr:9,assum:[1,2,3,4,5,13],assumd:4,astyp:[4,5,11],asum:4,asyp:11,atla:[0,2,7,12,13],atlas:1,author:[0,4,5,6,9,10],autmat:6,automat:[2,4,5,7,10],averag:6,avoid:5,axi:[1,5,7],background:[5,6,7],backgroundremov:0,backgrounds:[5,7],bain:2,ball:2,base:[2,4,5,6,9],basic:[5,11],becaus:[5,6],befor:7,belong:5,below:5,benjamini:3,best:[5,11],better:11,between:[2,7],bgr:11,bin:1,binari:[1,4,5,13],bonferoni:3,bool:[1,3,4,5,6,9,10],both:[4,13],boundari:5,box:5,brain:[0,2,5,11,12,13],brainactivitymap:[2,5,8,11],broad:5,bsplineparameterfil:[1,7],built:[4,5,7,10,11,12],button:5,calcul:[1,4,5,6,7,9],calculatechunks:5,calculatestatist:[5,11],calculatestatisticsonstack:[5,11],calculatestatisticsparamet:5,calculatesubstack:5,call:4,can:[1,2,3,4,5,6,11,12,13],care:6,cchunkoverlap:7,cell:[2,4,5,6,7],celldetect:0,cellintens:11,cellintensitymethod:5,cellprofil:5,cellshap:5,cellshape_d:5,cellshapedetect:11,cellshapethreshold:5,cellsiz:11,cellsizedetect:0,cellular:13,center:[2,4,5,6,11,13],certain:11,cfo:[5,11],chang:[1,2,5],channel:[4,13],chapter:5,check:[1,4],checkelastixiniti:1,chirstoph:[4,5,6],choos:[1,5],choosen:[5,11],christoph:[0,4,9,10,12],chunck:5,chunk:[5,7,11],chunkoptim:[5,7],chunkoptimizations:[5,7],chunkoverlap:[5,7],chunksiz:[5,11],chunksizemax:[5,7],chunksizemin:[5,7],cii:12,citi:[0,4,5,6,9,10],ckirst:[1,2,5,8,11],classic:13,classif:[5,7,13],classifi:[5,7],classifycel:5,classifycellsparamet:5,classifypixel:5,classifypixelparamet:5,cleamap:[4,13],cleanup:[1,2],clear:[0,1,2,5],clearanc:13,clearmappath:[1,4,5,8,11],clone:12,close:4,closefil:4,code:[4,6,11],col:2,collaps:2,collapsemap:2,collect:5,collumn:2,colon:4,color:[2,4,10,11,13],colorpvalu:2,column:[2,4],com:[3,12],combin:[2,5,6],comma:[4,13],comparison:[2,3],compat:13,complet:2,complex64:6,complex:13,compos:1,composit:4,comput:[1,3],concaten:4,configur:1,consid:5,consist:[1,4,5,13],consistentci:2,consit:9,constant:[2,3],constantli:3,consumpt:[5,6],contain:[1,11],continu:3,contribut:5,convert:[1,2,4,5,6],convertdata:4,convolut:5,convolv:6,coordiant:[1,4],coordin:[1,2,4,5,11,13],cooridn:2,coorind:1,copi:[2,4],copydata:4,copyfil:4,copyright:[4,5],cordin:11,corect:5,coron:1,correct:[2,3,5,6,13],correctillumin:5,correctilluminationparamet:[5,7],correctoin:7,correctpvalu:3,correspond:[1,4,5,13],cortic:13,could:5,count:2,countpointsgroupinregion:2,countpointsinregion:2,coupl:11,cramer:3,creat:[1,2,4,5,6],createdirectori:4,cretain:2,cross:2,csv:[0,2],csvfile:4,cubic:3,cumul:2,current:1,cutoffpvalu:2,cv2:[1,5],data:[0,1,2,3],databgr:11,datadog:11,datafileextens:4,datafileextensiontotyp:4,datafilenametomodul:4,datafilenametotyp:4,datafiletyp:4,datamax:11,datas:[1,2,4],dataset:2,datashap:11,datasizefromdatarang:4,datasizesink:1,datasizesourc:1,datasouc:10,datasourc:10,datatorang:4,datazs:4,date:4,david:4,dealt:4,deault:9,decreas:5,decrib:13,defaultannotationfil:2,defaultflatfieldlinefil:5,defaultlabeledimagefil:2,defin:[5,6,7,9],deform:1,deformationdist:1,deformationfield:1,delta:4,denot:4,densiti:7,describ:[3,4,5],descript:5,descrit:[2,4,5,6,13],descritpt:[2,5,6],desing:[4,5],destin:[1,4,5,10],detail:[4,12],detect:[2,5,6,7],detectcel:5,detectcellparamet:7,detectcellshap:[5,11],detectcellshapeparamet:[5,7],detectspot:5,detectspotparamet:5,detectspotsparamet:5,determin:[1,2,4,5,6,11],detremin:1,dict:[2,4,5,6,7,9],dictionari:[4,5,9],dictonari:9,differ:[2,3,4,5,6,11],differnet:5,digit:4,dilat:5,dim:5,dimens:[1,4,5],dimension:13,direct:[1,5],directli:[4,6],directori:[1,4,7],directorti:1,directri:1,discrimin:4,discuss:5,disk:[5,6,7],displai:1,distanc:1,distribut:[2,3,4,5],document:[1,12,13],doe:4,dog:[5,6],dogfilt:5,dogsiz:5,done:[1,5,13],dont:[5,6],dopaminerg:13,downhil:5,download:[8,12],drawn:3,dure:3,each:[1,2,4,5,9,11],earli:13,easier:11,easir:4,edit:[5,8],effect:5,effici:[4,5],either:[3,4,5],elaps:[5,9,11],elapsedtim:9,elast:1,elastix:0,elastixbinari:1,elastixlib:1,elastixparameteraffin:7,elastixparameterbsplin:7,elastixpath:[1,8],electron:5,element:5,elment:6,embryon:13,enabl:[4,13],endogen:13,entir:[4,5,11],equal_var:2,erod:5,eros:5,error:[4,5,13],espect:[2,13],essenc:1,essenti:5,estim:[1,3,5,11,13],estimateqvalu:3,etc:7,even:[6,11],eventhough:5,evert:4,evidenc:13,exampl:[1,4,5],except:4,exchang:1,exchng:7,exclud:2,execut:[1,3,5],exist:4,exmapl:[1,5],exp:4,expect:3,experi:11,experiment:4,express:[4,5],exten:5,extend:[3,4,5,6],extendedmax:5,extens:[4,5],extent:4,extern:0,exterpol:3,extra:[2,3,4],extract:[1,5],extrap1d:3,extrap1dfrominterp1d:3,extrapol:2,extrem:3,facil:13,facilit:13,factor:[1,5],fals:[2,3,4,5,6,7,9,10],fast:[5,13],fdr:3,featur:[3,5,11],felxibl:13,field:[1,2,4,5],figur:[10,11],file:[1,2,4,5,6,7,11,13],fileexperssiontofilenam:4,fileextens:4,filelist:0,filenam:[1,2,4,5,11],filter:5,filterdog:[5,6,11],filterdogparamet:[5,6,7],filterkernel2d:6,filterkernel3d:6,filterkernel:5,filterlinear:6,filterlinearparamet:6,find:[1,5],findcellintens:[5,11],findcellintensityparamet:5,findcells:[5,11],findcelsizeparamet:5,findcenterofmaxima:[5,11],findcenterofmaximaparamet:5,findextendedmaxima:[5,11],findextendedmaximaparamet:[5,7],findintens:5,findintensityparamet:[5,7],findpixelcoordin:5,first:[1,4,5,6,9,11],fit:[5,11],fitler:6,fix:[1,4],fixedimag:[1,7],fixedimagemask:7,fixinterpol:1,fixorient:1,flag:4,flat:5,flatfield:[5,7],flatfieldfromlin:5,flatfieldlinefromregress:5,flexibl:[4,5],float32:6,focu:13,folder:[1,4,5,8,13],follow:[4,5],form:[4,5,11,12],format:[1,2,4,5,9,13],formatelapsedtim:9,found:[1,2,5,13],four:[4,13],frame:1,framework:4,freeli:5,fro:4,from:[0,1,2,3,4,5,6,9,11,12],ftype:6,fuction:5,full:[4,5],fulli:5,func:5,functon:5,fundament:5,further:4,futur:5,fwer:3,gather:[5,11],gaussian:[5,6,7],gene:13,gener:[0,1,2,3,4,5,6,9,12],get:[4,5,9,11],getdataext:4,getdatas:4,getparamet:9,getresultdatafil:1,getscaleandoffset:4,gettransformfilesizeandspac:1,gettransformparameterfil:1,gfp:13,git:12,github:[3,12],give:[1,5],given:[1,2,3,4,5,6],gnu:[0,12],grai:[7,10],grayscal:5,grei:5,greyreconstruct:0,greyreconstructionparamet:5,greyscal:5,group1:2,group2:2,group:2,guassian:[6,11],h5file:4,h5py:4,hammond:4,handl:[4,5,9,10],have:[1,2],hdf5:4,head:9,header:[4,9,13],help:12,hemispher:7,here:[1,4,11,13],higher:5,hippocamp:13,hire:11,histolog:13,hmax:[5,7,11],hmaxtransform:5,hochberg:3,hold:2,holder:4,home:[1,2,5,8,11],horizont:11,host:8,hour:9,how:2,html:[3,4],http:[3,4,13],hypothesi:3,iamg:5,identifi:[1,5],idisco:[0,5],ieee:5,iid:2,ilasstik:5,ilast:7,ilastik:[5,7,8,13],ilastikclassif:0,ilastikparamet:7,ilastikpath:8,illastik:5,illumin:[5,13],illuminationcorrect:0,illuminationsc:5,illustr:11,imag:0,image:13,imageanalysi:11,imagej:[1,4,11,13],imageprocess:[0,4],imagestatist:0,imari:[0,2],img:[1,4,5,6],imglabel:5,imgmax:5,immedi:13,immunolabel:13,implement:[4,6],implment:4,includ:[5,6,8,13],increas:[3,7],independ:3,index:[4,12],indic:[1,4,5],indicat:2,individu:[2,4,5],infer:4,info:[1,2,3,4,5,6,9,13],infomr:2,inform:[1,2,5,6,11],infront:1,inherit:2,initi:[1,2,5,11],initializeelastix:1,initializeilastik:5,inner:6,inp:4,input:9,insert:4,inspect:[5,11],instal:[5,8,12],instens:5,institut:5,int16:5,int32:4,intact:13,intefac:5,integ:[4,5],integr:[4,5,6,13],intens:[2,4,5,11,13],intensities1:2,intensities2:2,intensitit:4,intensitygroup:2,intensityrow:2,intensti:5,inter:3,interact:13,interfac:[1,4,5,13],intern:[4,8,13],interp1d:3,interpol:[1,3],interpret:12,intrins:1,introduc:11,introduct:[],invari:5,invers:[1,10,11],inverseorient:1,inversli:1,invert:[1,7,10],ipdirect:3,isbn:5,isdatafil:4,isfil:4,isfileexpress:4,isi:13,isotrop:2,ispointfil:4,iter:4,jet:10,join:[1,4,5,9,11],joinparamet:9,joinpoint:5,joinstatist:5,kamentski:5,kannan:4,keep:5,kei:[2,5,6,9,13],kernel:[5,6,13],keyvaluepair:4,kind:3,kirst:[0,4,5,6,9,10,12],label:0,labelatcollaps:2,labelatlevel:2,labelcolormap:10,labeledimag:2,labelimag:4,labelinfo:2,labelpoint:2,labelrecord:2,labelsourc:10,labeltoacronym:2,labeltocolor:2,labeltonam:2,labl:[2,5,10],lablel:5,laplacian:6,larg:[0,1,4,5,11,12,13],larger:[5,11],layer:11,ld_library_path:1,learn:13,least:4,lee:5,length:[1,4],letter:5,level:[1,2,4,5],lib:1,librari:[0,1,4],licens:0,lie:11,light:5,lightsheet:[0,5,12],lightsheet_flatfield_correct:5,like:[4,13],line:[5,9],linear:[1,3,6,7],linearfilt:5,linearli:3,linux:1,list:[2,3,4,5,9,13],load:[4,11],local:[5,11,13],localmax:5,localmaxima:5,localmaxs:5,locat:[1,2],log:6,low:3,lowmemori:3,m_0:3,maarten:4,mai:[1,11],main:[1,4,5,11,13],make:[2,4,5],makecolorannot:2,makecolorpalett:2,mani:[4,5],map:[1,2,4,10,13],mark:2,mask:5,massachusett:5,match:[1,4],math:1,matlab:4,matplotlib:[1,11],max:[2,4,5,7,11],maxim:[2,5,7],maxima:5,maximadetect:0,maximalnumb:10,maximum:5,maxlabel:[5,11],maxtil:10,mean:[2,5,6,7,11],measur:[4,5,11],memori:[3,5,6,11],meta:4,meta_dict:4,method:[1,2,3,4,5,7],mhd:[1,4,13],microscop:[4,5],microscopi:[4,5,13],might:[2,5,11],miht:11,min:[4,5],minial:5,minim:[1,5,7],minimum:5,minu:1,minut:9,mise:3,mod:6,mode:[1,4,5,6],model:13,modifi:[3,4,5,6],modul:0,modular:[5,13],moment:4,more:[3,5,11,12,13],morpholog:[5,6,13],most:[5,11],mous:[0,2,5,11,12,13],move:[1,4],movingimag:[1,7],mpl:11,multi:4,multicomp:3,multipl:[1,3,4,9],multiplecomparisoncorrect:2,multipletest:3,must:[4,5],name:[1,2,4,5,6],ndarrai:3,nearest:1,nearli:[4,13],necessari:1,need:[1,2,4,5],neg:[1,2],negativetrend:2,net:[3,4],neuron:13,next:4,nfusi:3,nice:2,non:[1,5,7],none:[1,2,3,4,5,6,7,9,10,11],nonrigid:13,noprocess:5,note:[1,2,3,4,5,8,11],now:11,npy:0,nrdh:4,nrrd0005:4,nrrd:0,nrrderror:4,nrrdfile:4,nstack:5,nuclei:5,numbeer:1,number:[1,2,3,4,5,7,9,11,13],numer:4,numpi:[2,3,4,13],object:[1,2,3,4,5,6,9,10],observ:3,obtain:5,occur:9,odd:6,offset:[2,4,5,6],often:11,onc:[1,11],onli:[1,4,5,11],onlin:5,onot:5,onto:[1,13],open:[1,4,5,6,7,13],openfil:4,oper:[5,6],opject:4,optim:[4,5,6,7],optimizait:7,option:[2,3,4,5,11,13],order:[1,3,4,5],ordereddict:2,org:[12,13],organ:[9,13],orientationtopermu:1,orientdatas:1,orientdatasizeinvers:1,orientresolut:1,orientresolutioninvers:1,origin:[1,4,5,6],ostenrefara_v2_lowerhalf:7,other:[4,5,7,13],otherwis:[1,5],our:11,out:[1,2,5,6,9],outer:6,output:[1,2,4,7,9],outsid:3,overlai:[10,11],overlap:[5,7],overlaylabel:10,overlaypoint:10,overrid:4,own:5,page:12,pair:[4,9],pal:2,paral:5,parallel:[0,1,4],parallelprocessstack:5,paralllel:5,paramet:[0,1,2,3,4,5,6],parametertool:0,paramt:7,parent:2,pars:[1,4],parse_nrrdvector:4,parse_optional_nrrdvector:4,parseelastixoutputpoint:1,part:[2,5,13],particular:[1,2,3,11],pass:4,path:[1,4,5,8,11],pattern:[4,5,13],patttern:4,pcutoff:2,pdf:11,peak:[5,13],per:1,percentil:2,perform:[1,2,5,6,11],permuat:1,permut:1,physic:1,pi0:3,pickl:2,pipelin:5,piplin:5,pixel:[1,2,4,5,7],place:[4,5],placehold:4,plain:2,plane:11,platform:4,plot:[0,2,5],plotoverlaylabel:[10,11],plotoverlaypoint:[10,11],plottil:[10,11],plt:11,pmax:2,png:11,point:[0,1,2],pointcolor:10,pointcounts1:2,pointcounts2:2,pointfileextens:4,pointfileextensiontotyp:4,pointfilenametomodul:4,pointfilenametotyp:4,pointfiletyp:4,pointgroup:2,points1:2,points2:[2,4],pointshiftfromrang:4,pointsink:1,pointsourc:[1,10],pointstocoordin:4,pointstocoordinatesandproperti:4,pointstocoordinatesandpropertiesfilenam:4,pointstoproperti:4,pointstorang:4,polat:3,polynomi:5,porcess:5,posit:[1,2,3,5],positivetrend:2,possibl:[4,11,13],postfix:4,potenti:[1,5],predefin:5,prefix:9,premut:1,prep:13,present:6,press:5,primari:1,principl:5,print:[1,3,4,5,6,9,11],printelapsedtim:9,printset:1,printsubstackinfo:5,probabl:5,proces:5,process:[0,1,4],processig:5,processingdirectori:1,processmethod:5,processstackparamet:7,processwrit:0,profil:5,program:[1,8],progress:[1,5,6,11],project:[2,5,8,11],propabl:5,properli:5,properti:4,propertiespostfix:4,proport:3,protocol:4,provid:[1,2,4,5,6,9,10,11,13],psign:2,pull:3,pval:2,pvalu:3,pyplot:11,python:[0,3,4,12],quantifi:13,question:2,qvalu:3,radial:5,radiu:[2,4,6],rand:4,random:4,rang:[2,4,5,10,11,13],rapid:13,raster:[4,13],ratio:3,raw:[0,1],read:[4,11,13],readdata:[4,5,11],readdatafil:4,readdatagroup:2,readdataset:4,reader:4,readfilelist:4,readhead:4,readpoint:4,readpointsgroup:2,real:4,recognit:5,reconstruct:[5,13],rectangular:[2,7],redesign:5,reduc:4,reduct:4,redund:5,refer:[1,2,3,5,7],referenec:[2,13],refin:5,reg:4,region:[2,5,13],regist:[1,13],registr:[0,1,12,13],regoin:2,regular:4,rejoin:5,rel:[1,4],remov:[1,5],removebackground:[5,11],removebackgroundparamet:[5,7],removenan:2,render:4,renier2014:13,renier2015:13,renier:13,reorient:[1,13],replac:1,report:13,repres:[4,13],represen:1,requir:5,resampl:0,resampledata:1,resampledatainvers:1,resampledatas:1,resamplepoint:1,resamplepointsinvers:1,resamplexi:1,resamplingparamet:7,rescal:[1,5,7],rescalefactorilastik:5,rescalesizeandspac:1,rescaletoilastik:5,reset:9,reslut:2,resmapl:1,resolut:[1,4,7,13],resolutionsink:[1,7],resolutionsourc:[1,7],respect:[2,13],restrict:4,result:[1,2,4,5,6,7,9,10,11],resultdir:1,resultdirectori:1,returncount:2,returnid:2,revers:[1,2,5],rewrit:4,rgb:10,rigid:13,robinson:5,rockefel:[0,4,5,6,9,10,12],root:[1,8],routin:[1,4,5,6,8,10,11],row:2,run:[1,5,11],saggit:1,sagittaltocoronaldata:1,same:[2,3,6],sampl:[0,1,3,4,5,12,13],sandbox:3,save:[1,4,5,6,7],scale:[1,4,5,7,10],scienc:[2,5,8,11],scientif:2,scikit:[3,5],scipi:[3,5,6],script:5,search:[5,12],second:[4,6,9],secondari:1,section:[10,12,13],see:[3,4,5,6,11,12,13],seed:5,segment:[11,13],selem:5,self:2,separ:[4,13],separatehead:4,sequenc:[1,3],sequenti:5,sequentiallyprocessstack:5,sequnec:10,sesiz:6,set:[0,1,2,4,5,6,7],setelastixlibrarypath:1,setpathtransformparameterfil:1,sett:8,settransformfilesizeandspac:1,setup:[1,8],setyp:6,sever:[4,5,13],shape:[1,4,5],sheet:5,shift:[4,5],shiftpoint:5,should:[1,5],shown:11,side:3,sigma2:[6,7],sigma:[6,7],sign:[1,2],signal:6,signific:2,simpl:[4,9,13],simpli:4,sinc:9,singl:[4,5],sink:[1,2,4,5,7,10],slf:2,slice:[1,4,5],small:11,smaller:[2,5],smear:2,softwar:[8,11,13],soill:5,solvent:13,some:[2,3,5,11],sometim:11,sort:2,sought:1,sourc:[1,2,3,4,5,6,7,8,9,10,11],sourceforg:[3,4],space:[1,4],spatial:[1,4],special:[4,6,8,11],specif:[1,2,4,5,6,8,9,10,13],specifi:[1,4,5,11],sped:6,speed:[5,6],sphere:[4,6],spheric:[2,7],spline:3,split:[4,5],splitfileexpress:4,spot:[4,5,7],spotdetect:0,spotdetectionparamet:7,stack:[1,2,4],stackid:5,stackprocess:0,stackprocessingparamet:7,stain:[5,11,13],stand:13,standard:[1,3,4,5],start:[5,9],startindex:4,stat:[3,11],statist:0,statisticaltest:2,statisticsfrom:5,statsmodel:3,std:[2,6],stdout:[1,5,6,9],step:[1,5,6,11],still:4,stop:9,store:[1,4,13],storei:3,str:[1,2,3,4,5,6,8,9,10],string:[1,2,4,9],structur:[2,5],structureel:5,structureelement2d:6,structureelement3d:6,structureelementoffset:6,sub:[1,2,4],subregion:11,subsequ:1,subset:4,substack:[5,6,11],subtract:5,success:4,sucessfulli:5,suitabl:4,sum:5,summari:0,support:[1,2,4,5,10,13],sure:4,synthet:5,system:[5,13],systemat:5,tabl:[1,2,4],tail:3,take:[1,5,6],taken:5,target:[0,1,12,13],techniqu:13,technolog:5,teem:4,temporari:1,test:[1,2,3,4,5,6,7,10,11],test_idisco_:5,testcompletedcumul:2,testcompletedcumulativesinspher:2,testcompletedinvertedcumul:2,testcramervonmises2sampl:3,text:[4,9,13],than:[3,5],thei:4,them:11,therfor:11,thi:[1,2,3,4,5,6,7,8,9,10,11,13],third:4,those:[5,11],three:13,threshold:[2,5,7,11],thresholdpoint:2,throughout:9,thu:5,tibshirani:3,tif:[0,1,2],tiff:4,tifffil:4,tile:10,time:[4,5,9,11],timepoint:4,timer:0,tissu:[0,12,13],tmpfile:1,todatarang:4,todatas:4,todo:[2,4],tolabelatcollaps:2,tolabelatcollapsemap:2,tolabelatlevel:2,tomultichanneldata:4,too:5,tool:[1,2],toolbox:[0,2,4,5,12,13],toolkit:13,top:1,total:[3,5,11],toward:[0,12,13],train:5,trale:4,trane:5,transact:5,transform:[1,4,5],transformdata:1,transformdirectori:1,transformfil:1,transformix:1,transformixbinari:1,transformparameterfil:1,transformpoint:1,transformpointstoimari:4,transfrom:5,transgen:13,translat:13,transpar:10,trivial:5,ttestpointsinregion:2,ttestvoxel:2,tupl:[1,2,4,5,6,10],tuplr:1,turn:2,two:[3,4,5,7],txt:[4,7],type:[1,2,4,5],typic:[1,5],uin8:5,uint8:5,umadevi:4,ummber:6,uniform:[2,6],univers:[0,4,5,6,9,10,12],usag:11,useful:[4,11],user:5,usual:1,util:[0,4,6],valid:4,valu:[2,3,4,5,9,13],variabl:[1,4],variat:[5,11],variou:[2,3,4,6,7],vebos:1,vector:[1,4,5],venkataraju:4,verbos:[1,3,5,6,7,9,11],veri:[4,5],version:[0,3,5,12],via:[1,3,5,10,11,13],video:13,view:[1,4,11,13],viewer:[4,13],vignet:[5,13],vincent:5,visibl:11,visual:[0,2,4],visula:2,volum:[5,7,13],volumetr:[0,1,2,4],von:3,voxel:0,voxelizationparamet:7,voxelizations:7,voxelizeparamet:2,voxelizepixel:2,vtk:0,wai:[1,5,9,13],warpabl:7,watersh:[5,11,13],webpag:[8,13],weight:[2,5],weightsfromprecentil:2,were:[2,5],when:[1,4,5,11,13],where:[1,4,5,13],which:[1,2,3,4,5,7],wide:[2,13],window:4,without:[1,2,4],work:[1,4,5],would:4,write:[1,2,4,5,6,9,13],writedata:4,writehead:4,writelut:2,writep:2,writeparamet:9,writepoint:[1,4],writer:4,writerawdata:4,writestr:9,writesubstack:5,writet:4,written:4,wrt:7,www:13,xlabel:11,xsize:5,yet:4,ylabel:11,york:[0,4,5,6,9,10,12],you:12,zcenter:5,zcenterindic:5,zero:[2,5],zoom:11,zsubstackcenterindic:5},titles:["ClearMap package","ClearMap.Alignment package","ClearMap.Analysis package","ClearMap.Analysis.Tools package","ClearMap.IO package","ClearMap.ImageProcessing package","ClearMap.ImageProcessing.Filter package","ClearMap.Parameter module","ClearMap.Settings module","ClearMap.Utils package","ClearMap.Visualization package","ClearMap Image Analysis Tools","ClearMap","Introduction"],titleterms:{align:[1,13],analysi:[2,3,11,13],author:12,background:11,backgroundremov:5,cell:11,celldetect:5,cellsizedetect:5,clear:13,clearmap:[0,1,2,3,4,5,6,7,8,9,10,11,12],content:12,convolut:6,csv:4,data:4,detect:11,dogfilt:6,elastix:1,element:6,exampl:13,extern:5,extrapol:3,filelist:4,filter:[6,11],filterkernel:6,greyreconstruct:5,idisco:13,ilastikclassif:5,illuminationcorrect:5,imag:[1,4,5,11,13],imageprocess:[5,6],imagestatist:5,imari:4,indic:12,introduct:13,label:2,licens:12,linearfilt:6,maxima:11,maximadetect:5,method:13,modul:[1,2,3,4,5,6,7,8,9,10],multiplecomparisoncorrect:3,npy:4,nrrd:4,orient:1,packag:[0,1,2,3,4,5,6,9,10],parallel:5,paramet:7,parametertool:9,plot:[10,11],point:4,process:[5,13],processwrit:9,quickstart:12,raw:4,refer:13,remov:11,represent:1,resampl:1,set:8,shape:11,size:1,spotdetect:5,stack:5,stackprocess:5,statist:[2,11],statisticaltest:3,structur:6,structureel:6,sub:5,subpackag:[0,2,5],summari:[1,4],tif:4,tile:11,timer:9,tool:[3,11],type:6,util:9,visual:[10,11],volumetr:5,voxel:2,vtk:4}})