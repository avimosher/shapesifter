import os

#options=Options('SConstruct.options')

external_libraries_dir="#External_Libraries/"
external_libraries={
    'eigen': {'default': 1, 'libs':[''],'cpppath':[external_libraries_dir+'eigen',external_libraries_dir+'eigen/unsupported']},
    'json': {'default': 1,'cpppath':[external_libraries_dir+'jsoncpp/dist'],'libs':['jsoncpp'],'libpath':[external_libraries_dir+'jsoncpp/dist']}
}

env=Environment()

base_env=Environment()

builder=base_env.SharedObject
jsoncpp='#External_Libraries/jsoncpp/dist/jsoncpp.cpp'
jsoncpp_obj=builder(target=os.path.splitext(jsoncpp)[0],source=jsoncpp)
base_env.SharedLibrary(target='#External_Libraries/jsoncpp/dist/jsoncpp',source=jsoncpp_obj)


for name,lib in external_libraries.items():
    defaults={'default':0,'flags':'','linkflags':'','cpppath':[],'libpath':[]}
    for f in defaults.keys(): lib.setdefault(f,defaults[f])
    env.Append(CPPPATH=lib['cpppath'])
    env.Append(LIBPATH=lib['libpath'])
    env.Append(LIBS=lib['libs'])
    env.Append(LINKFLAGS=lib['linkflags'])
#    options.AddOptions(BoolOption('USE_'+name.upper(),'Use '+name,lib['default']),
#                       (name+'_include','Include directory for '+name,0),
#                       (name+'_libpath','Library directory for '+name,0),
#                       (name+'_rpath','Extra rpath directory for '+name,0),
#                       (name+'_libs','Libraries for '+name,0),
#                       (name+'_linkflags','Linker flags for '+name,0))

#conf=Configure(env)
#print(conf.CheckLib('jsoncpp'))

def Automatic_Program(target,source,env):
    env.Program(target=target,source=source)
    env.Install('#bin',target)

def Automatic_Library(target,source,env):
    library=env.SharedLibrary(target=target,source=source)
    print(library)
    env.Install('#bin',library)

env.Append(CCFLAGS="-std=c++11 -g")
env.Append(CPPPATH="#Library")
directories=SConscript('Library/SConscript',variant_dir='build/Library',exports={'env': env,'Automatic_Library': Automatic_Library})

env_projects=env.Copy()
env_projects.Append(LIBPATH=['#build/Library'])
env_projects.Append(LIBS=directories)
#print(libraries)
#env.Install('#bin',libraries)

SConscript('Projects/SConscript',variant_dir='build/Projects',
                       exports={'env': env_projects,
                                'Automatic_Program': Automatic_Program})
#print(executables)
#env_projects.Install('#bin',executables)
#SConscript('Tests/SConscript',variant_dir='build/Tests')
