# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 19:51:24 2015

@author: mtllab
"""


#basepar = '/home/mtllab/Documents/warping/Par0000bspline_test_';

parfiles = ['5', '5_a', '5_b', '5_c', '5_d', '5_e', '5_f',
            '25', '25_a', '25_b', '25_c', '25_d', '25_e', '25_f',
            'a', 'a_a', 'a_b', 'a_c'];
                  
execfile('1266.py');
        
for p in parfiles:
    
    RegistrationAlignmentParameter["resultDirectory"] = os.path.join(BaseDirectory, 'elastix_auto_to_atlas_test_' + p);
    RegistrationAlignmentParameter["bSplineParameterFile"] = os.path.join(PathReg, 'Par0000bspline_test_' + p + '.txt');
    
    resultDirectory  = alignData(**RegistrationAlignmentParameter);



execfile('dulac1.py');
        
for p in parfiles:
    
    RegistrationAlignmentParameter["resultDirectory"] = os.path.join(BaseDirectory, 'elastix_auto_to_atlas_test_' + p);
    RegistrationAlignmentParameter["bSplineParameterFile"] = os.path.join(PathReg, 'Par0000bspline_test_' + p + '.txt');
    
    resultDirectory  = alignData(**RegistrationAlignmentParameter);
