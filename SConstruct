import glob, os

env = Environment()

include_path = [#'/home/kevin/armadillo/include',
				#'/home/kevin/pca/include',
				'#.' ]

env.Append( CPPPATH = include_path )

lib_path = [#'/home/kevin/armadillo/lib',
			#'/home/kevin/pca/lib'
			]

env.Append(LIBPATH = lib_path)


libs = [#'pca',
		'boost_system',
		'tbb',
		'flann_cpp_s',
		'OpenImageIO',
		'Imath']

env.Append( CCFLAGS = '-O0 -g -fopenmp -std=c++0x' )
env.Append( LIBS = libs)

cpp_defines = [ ]
env.Append( CPPDEFINES = cpp_defines)

env.Append( LINKFLAGS = '-fopenmp' )

sourceFiles = Glob( '*.cpp' )
sourceFiles += Glob( './analyzer/*.cpp' )
sourceFiles += Glob( './synthesizer/*.cpp' )

env.Program('texsyn', sourceFiles)

