#!/usr/bin/python

from optparse import OptionParser
from jinja2 import Template
from subprocess import call
from tempfile import mkstemp
import os
import ast
import itertools
from datetime import date, timedelta, datetime
import glob
import csv
import numpy
import json
from lxml import etree
import math

parser = OptionParser()
#parser.add_option("-f",dest="filename",type="string")
parser.add_option("-m",action="store_true",dest="mean",default="False")
(options, args) = parser.parse_args()

output_types=1
results_structure={}
results_structure['results']=[]
expected_samples=0

def analyze(path,prefix,parameter_values,fraction_done):
    global output_types
    global observed_samples
    global results_structure

    pattern=path+"/"+prefix+".[0-9]*"
    
    count=0
    total_quantity=0
    values=[]
    for name in glob.glob(pattern):
        with open(name, 'r') as csvfile:
            csvreader = csv.reader(csvfile, delimiter='\t')
            for row in csvreader:
                try:
                    value_tuple=[]
                    for v in row:
                        try:
                            value_tuple.append(float(v))
                        except:
                            pass
                    values.append(value_tuple)
                    count=count+1
                except:
                    pass
    result={}
    result['parameters']=parameter_values
    result['expected']=expected_samples
    result['observed']=count
    if options.mean!=True:
        result['values']=values
    result['fractionDone']=fraction_done
    meanValues=tuple(i/count for i in [sum(x) for x in zip(*values)])
    deviation=[]
    for i in range(0,len(meanValues)):
        total_deviation=0
        for j in range(0,len(values)):
            value=values[j]
            total_deviation+=math.pow(meanValues[i]-value[i],2)
        deviation.append(math.sqrt(total_deviation)/count)
    result['deviation']=deviation
    result['meanValues']=meanValues
    results_structure['results'].append(result)
    if count>0:
        output_types=len(values[0])
    return values



def add_results(filename):
    global expected_samples
    global results_structure
    parameter_file=open(filename,'r')

    ast_string=parameter_file.read()

    control_object=ast.literal_eval(ast_string)

    scene_template_filename=control_object['scene']
    f=open(scene_template_filename, 'r')
    scene_template_string=f.read()

    parameters_list=[]
    key_list=[]
    for key in control_object['parameters']:
        key_list.append(key)
        parameters_list.append(control_object['parameters'][key])

    column_label=control_object['column_label']
    expected_samples=control_object['iterations']
    results_structure['variables']=key_list
    results_structure['variableValues']=parameters_list

    for a in itertools.product(*parameters_list):
        parameter_dict={}
        name_list=[]
        parameter_values=[]
        for i in range(0,len(key_list)):
            parameter_dict[key_list[i]]=a[i]
            name_list.append(key_list[i])
            name_list.append(str(a[i]))
            parameter_values.append(float(a[i]))
            
        template=Template(scene_template_string)
        scene_string=template.render(parameter_dict)
        root=etree.fromstring('<data>\n'+scene_string+'\n</data>')

        #last_time=float(root.find('example').find('last_time').text)
        #minimum_dt=float(root.find('example').find('minimum_dt').text)
        #output_directory=root.find('example').find('output_directory').text
        #log_file=output_directory+'/common/log.txt'
        disable_candidate=False
        #try:
        #    logroot=etree.parse(log_file)
        #except etree.XMLSyntaxError:
        #    disable_candidate=True
        #maximum_frame=last_time/minimum_dt
        #current_frame=0
        ## count down from end
        ## nasty hack, fix by changing output from simulation instead
        #for i in range(1,20):
        #    if disable_candidate==False:
        #        candidate_frame=logroot.xpath('./print[last()-'+str(i)+']')[0].text
        #        elements=candidate_frame.split(' ')
        #        if len(elements)>=4 and elements[0]=='Current':
        #            current_frame=int(elements[3])
        #            break
        #    else:
        #        current_frame=maximum_frame
        output_directory="output/"+control_object['prefix']+'_'+'_'.join(name_list)
        #print output_directory
        analyze(output_directory,"output",tuple(parameter_values),1)#current_frame/maximum_frame)

for filename in args:
    add_results(filename)

results_structure['retrievedTime']=datetime.now().isoformat()
print json.dumps(results_structure,sort_keys=False,indent=4,separators=(',',': '))

