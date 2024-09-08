%-------------------
% Compile Mex files.
%-------------------
% cd multigs;
% if exist('computeIntersection','file')~=3
%     mex computeIntersection.c;
% end
% cd ..;
% %% mex file for warping 2017.3.30
% if exist('imageWarping','file')~=2 || exist('imageProjection','file')~=2  || exist('getLocation','file')~=2 ||...
%    exist('getBlendMask','file')~=2 || exist('ceresRigidError','file')~=2 || exist('ceresNonrigidError','file')~=2
%     mex imageWarping.cpp;
%     mex imageProjection.cpp;
%     mex getBlendMask.cpp;
%     mex getLocation.cpp;
%     
%    % CERES solver and GLOG (both from google) should be installed in your
%    % system in order to compile the following files.
%    % mex ceresRigidError.cpp /usr/local/lib/libceres_shared.so /usr/local/lib/libglog.so -I/usr/include/eigen3;
%    % mex ceresNonrigidError.cpp /usr/local/lib/libceres_shared.so /usr/local/lib/libglog.so -I/usr/include/eigen3;
%    mex LINKFLAGS="$LINKFLAGS /LIBPATH:D:\3rd\x64\Release ceres.lib libglog_static.lib" ...
%    COMPFLAGS="$COMPFLAGS /openmp ...
%    ... /D DEBUG ...
%    /D GOOGLE_GLOG_DLL_DECL= ...
%    /I D:\3rd\ceres-solver\include /I D:\3rd\win\include ...
%    /I D:\3rd\google-glog\src /I D:\3rd\google-glog\src\windows /I D:\3rd\Eigen" ...
%    ceresRigidError.cpp
%    % copyfile('D:\3rd\x64\Release\ceres.dll','.')
%    mex LINKFLAGS="$LINKFLAGS /LIBPATH:D:\3rd\x64\Release ceres.lib libglog_static.lib" ...
%    COMPFLAGS="$COMPFLAGS /openmp ...
%    ... /D DEBUG ...
%    /D GOOGLE_GLOG_DLL_DECL= ...
%    /I D:\3rd\ceres-solver\include /I D:\3rd\win\include ...
%    /I D:\3rd\google-glog\src /I D:\3rd\google-glog\src\windows /I D:\3rd\Eigen" ...
%    ceresNonrigidError.cpp
% ������
%    mex LINKFLAGS="$LINKFLAGS /LIBPATH:D:\Computer_Vision\3rd\x64\Release ceres.lib libglog_static.lib" ...
%    COMPFLAGS="$COMPFLAGS /openmp ...
%    ... /D DEBUG ...
%    /D GOOGLE_GLOG_DLL_DECL= ...
%    /I D:\Computer_Vision\3rd\ceres-solver\include /I D:\Computer_Vision\3rd\win\include ...
%    /I D:\Computer_Vision\3rd\google-glog\src /I D:\Computer_Vision\3rd\google-glog\src\windows /I D:\Computer_Vision\3rd\Eigen" ...
%    ceresRigidError.cpp

%    mex LINKFLAGS="$LINKFLAGS /LIBPATH:D:\Computer_Vision\3rd\x64\Release ceres.lib libglog_static.lib" ...
%    COMPFLAGS="$COMPFLAGS /openmp ...
%    ... /D DEBUG ...
%    /D GOOGLE_GLOG_DLL_DECL= ...
%    /I D:\Computer_Vision\3rd\ceres-solver\include /I D:\Computer_Vision\3rd\win\include ...
%    /I D:\Computer_Vision\3rd\google-glog\src /I D:\Computer_Vision\3rd\google-glog\src\windows /I D:\Computer_Vision\3rd\Eigen" ...
%    ceresPointLineofLength_v1.cpp  


   mex LINKFLAGS="$LINKFLAGS /LIBPATH:D:\Computer_Vision\3rd\x64\Release ceres.lib libglog_static.lib" ...
   COMPFLAGS="$COMPFLAGS /openmp ...
   ... /D DEBUG ...
   /D GOOGLE_GLOG_DLL_DECL= ...
   /I D:\Computer_Vision\3rd\ceres-solver\include /I D:\Computer_Vision\3rd\win\include ...
   /I D:\Computer_Vision\3rd\google-glog\src /I D:\Computer_Vision\3rd\google-glog\src\windows /I D:\Computer_Vision\3rd\Eigen" ...
   ceresPointLineofLengthRef.cpp

