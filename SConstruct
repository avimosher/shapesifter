import os
import glob

Decider('timestamp-match')

variables=Variables(None,ARGUMENTS)
variables.Add(EnumVariable('TYPE','Type of build','release',allowed_values=('release','debug')))

external_libraries_dir="#External_Libraries/"
external_libraries={
    'catch': {'default': 1, 'libs':[''], 'cpppath':[external_libraries_dir+'catch']},
    'cereal': {'default': 1, 'libs':[''],'cpppath':[external_libraries_dir+'cereal/include']},
    'eigen': {'default': 1, 'libs':[''],'cpppath':[external_libraries_dir+'eigen',external_libraries_dir+'eigen/unsupported']},
    'ESBTL': {'default': 0, 'flags': ['NDEBUG'],'libs':['CGAL','mpfr','gmp'],'cpppath':[external_libraries_dir+'ESBTL/include']},
    'gmm': {'default': 1, 'libs':[''], 'cpppath':[external_libraries_dir+'gmm']},
    'json': {'default': 1,'cpppath':[external_libraries_dir+'jsoncpp/dist'],'libs':['jsoncpp'],'libpath':[external_libraries_dir+'jsoncpp/dist']},
    'osg': {'default': 0,'cpppath':[external_libraries_dir+'osg/include'],'libs':['osg','osgDB','osgGA','osgViewer','libOpenThreads','libosgUtil','libosgText'],'libpath':[external_libraries_dir+'osg/lib']},
    'gl': {'default': 0,'flags': ['VIEWER'],'libs':['GL']}
}

env=Environment(variables=variables)#,ENV={'PATH' : os.environ['PATH'], 'LD_LIBRARY_PATH' : os.environ['LD_LIBRARY_PATH']})
base_env=Environment()


if env['PLATFORM']=='darwin':
    env['PLATFORM_ALIAS']='OSX'
elif env['PLATFORM']=='posix':
    env['PLATFORM_ALIAS']='Linux'

build_base='build/'+env['PLATFORM']+'/'+env['TYPE']

base_env.Append(CPPPATH=external_libraries['json']['cpppath'])
builder=base_env.SharedObject
jsoncpp='#External_Libraries/jsoncpp/dist/jsoncpp.cpp'
jsoncpp_obj=builder(target=os.path.splitext(jsoncpp)[0],source=jsoncpp)
jsonlib=base_env.SharedLibrary(target='#External_Libraries/jsoncpp/dist/jsoncpp',source=jsoncpp_obj)
env.Install("#Shapesifter/"+env['PLATFORM_ALIAS'],jsonlib)
env.Install("#bin/"+env['PLATFORM_ALIAS'],jsonlib)


def Load_External(env):
    for name,lib in external_libraries.items():
        defaults={'default':0,'flags':'','linkflags':'','cpppath':[],'libpath':[],'flags':[]}
        for f in defaults.keys(): lib.setdefault(f,defaults[f])
        if env.get('USE_'+name.upper()) or lib['default']:
            env.Append(CPPDEFINES=lib['flags'])
            env.Append(CPPPATH=lib['cpppath'])
            env.Append(LIBPATH=lib['libpath'])
            env.Append(LIBS=lib['libs'])
            env.Append(LINKFLAGS=lib['linkflags'])

def Automatic_Program(target,source,env):
    local_env=env.Clone()
    Load_External(local_env)
    program=local_env.Program(target=target,source=source)
    local_env.Install('#bin/'+env['PLATFORM_ALIAS'],target)
    return target
    
def Automatic_Library(target,source,env):
    local_env=env.Clone()
    Load_External(local_env)
    library=local_env.SharedLibrary(target=target,source=source)
    local_env.Install('#bin/'+env['PLATFORM_ALIAS'],library)
    return library

def Find_SConscripts(env,dir):
    for c in glob.glob(os.path.join(dir,"*","SConscript")):
        env.SConscript(c,variant_dir=build_base+'/'+os.path.dirname(c),exports={'env': env,'Automatic_Program': Automatic_Program})
    for c in glob.glob(os.path.join(dir,"*","SConscript_unbuilt")):
        env.SConscript(c,exports={'env': env})

# Compiler flags
env.Append(CXXFLAGS="-std=c++1y -frounding-math")
if env['TYPE']=='debug':        
    env.Append(CXXFLAGS="-g")
elif env['TYPE']=='release':
    env.Append(CXXFLAGS="-O3")

env.Append(CPPPATH="#Library")
#env.Append(CXXFLAGS="-Wall -Winit-self -Woverloaded-virtual -Wstrict-aliasing=2 -fno-strict-aliasing -Wno-unused-but-set-variable -Werror")
directories=SConscript('Library/SConscript',variant_dir=build_base+'/Library',exports={'env': env,'Automatic_Library': Automatic_Library,'Load_External': Load_External})

env_projects=env.Clone()
env_projects.Append(LIBS=directories)
env_projects.Append(LIBPATH=['#'+build_base+'/Library'])

Find_SConscripts(env_projects,'Projects')
Find_SConscripts(env_projects,'Tests')
env.SConscript('bin/SConscript',exports={'env': env})


#############################################################
# Build release structure
#############################################################
