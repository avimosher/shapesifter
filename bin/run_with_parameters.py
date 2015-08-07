#!/usr/bin/python
from optparse import OptionParser
from jinja2 import Template
from subprocess import call
from tempfile import mkstemp
import os
import ast
import itertools

parser = OptionParser()
parser.add_option("-f",dest="filename",type="string")
parser.add_option("-t",dest="test",action="store_true")
(options, args) = parser.parse_args()

# read in control parameters
parameter_file=open(options.filename,'r')
ast_string=parameter_file.read()
print ast_string
control_object=ast.literal_eval(ast_string)

# get scene template file as string
scene_template_filename=control_object['scene']
f=open(scene_template_filename, 'r')
scene_template_string=f.read()

lo=open('log.txt','w')

parameters_list=[]
key_list=[]
for key in control_object['parameters']:
    key_list.append(key)
    parameters_list.append(control_object['parameters'][key])

def bsub_multiple(output_directory,iterations,prefix,output_name,queue,length,log_directory,parameter_dict):
    print output_directory
    template=Template(scene_template_string)
        
    for i in range(1,iterations+1):
        target_string=output_directory + "/" + output_name + "." + str(i)
        target_directory=output_directory + "." + str(i)
        command_target_string = prefix + "/" + target_string
        parameter_dict['output_directory']=log_directory+"/"+target_directory
        scene_string=template.render(parameter_dict)
        (f_out,scene_file)=mkstemp(dir='temporary_scripts')
        os.write(f_out,scene_string)

        command=control_object['command']%(scene_file)
        try:
            os.mkdir(log_directory+"/"+target_directory)
        except OSError:
            pass
        try:
            os.mkdir(prefix+"/"+output_directory)
        except OSError:
            pass

        full_command = "bsub -u null -q %s -W %s %s \> %s"%(queue,length,command,command_target_string)
        print full_command
        if not options.test:
            call(full_command, shell=True)

for a in itertools.product(*parameters_list):
    parameter_dict={}
    name_list=[]
    for i in range(0,len(key_list)):
        parameter_dict[key_list[i]]=a[i]
        name_list.append(key_list[i])
        name_list.append(str(a[i]))

    output_directory=control_object['prefix']+'_'+'_'.join(name_list)

    lo.write(os.getcwd())
    bsub_multiple(output_directory=output_directory,prefix="output",output_name="output",iterations=control_object['iterations'],queue=control_object['queue'],length=control_object['W'],log_directory="log_directory",parameter_dict=parameter_dict)

