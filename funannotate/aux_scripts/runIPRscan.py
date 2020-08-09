#!/usr/bin/env python
# $Id: iprscan5_urllib2.py 2809 2015-03-13 16:10:25Z uludag $
# ======================================================================
#
# Copyright 2009-2014 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# ======================================================================
# InterProScan 5 (REST) Python client using urllib2 and
# xmltramp (http://www.aaronsw.com/2002/xmltramp/).
#
# Tested with:
#  Python 2.6.5 (Ubuntu 10.04 LTS)
#  Python 2.7.3 (Ubuntu 12.04 LTS)
#
# See:
# http://www.ebi.ac.uk/Tools/webservices/services/pfa/iprscan5_rest
# http://www.ebi.ac.uk/Tools/webservices/tutorials/python
# ======================================================================
# Base URL for service
import urllib.request, urllib.error, urllib.parse
import urllib.request, urllib.parse, urllib.error
import time
import sys
import re
import os
import platform
import argparse
import xmltramp
baseUrl = 'http://www.ebi.ac.uk/Tools/services/rest/iprscan5'

# Load libraries

# Set interval for checking status
checkInterval = 10
# Output level
outputLevel = 1
# Debug level
debugLevel = 0
# Number of option arguments.
numOpts = len(sys.argv)

# Usage message
parser = argparse.ArgumentParser()
# Tool specific options
parser.add_argument('--input', required=True, help='input FASTA file')
parser.add_argument('--appl', 
					help='signature methods to use, see --paramDetail appl')
parser.add_argument('--crc', action="store_true",
                    help='enable InterProScan Matches look-up (ignored)')
parser.add_argument('--nocrc', action="store_true",
                    help='disable InterProScan Matches look-up (ignored)')
parser.add_argument('--goterms', action="store_true",
                    help='enable inclusion of GO terms')
parser.add_argument('--nogoterms', action="store_true",
                    help='disable inclusion of GO terms')
parser.add_argument('--pathways', action="store_true",
                    help='enable inclusion of pathway terms')
parser.add_argument('--nopathways', action="store_true",
                    help='disable inclusion of pathway terms')
parser.add_argument('--sequence', help='input sequence file name')
# General options
parser.add_argument('--email', required=True, help='e-mail address')
parser.add_argument('--title', help='job title')
parser.add_argument('--outfile', help='file name for results')
parser.add_argument('--outformat', help='output format for results')
parser.add_argument('--async', action='store_true', help='asynchronous mode')
parser.add_argument('--jobid', help='job identifier')
parser.add_argument('--polljob', action="store_true", help='get job result')
parser.add_argument('--status', action="store_true", help='get job status')
parser.add_argument('--resultTypes', action='store_true',
                    help='get result types')
parser.add_argument('--params', action='store_true',
                    help='list input parameters')
parser.add_argument('--paramDetail', help='get details for parameter')
parser.add_argument('--quiet', action='store_true',
                    help='decrease output level')
parser.add_argument('--verbose', action='store_true',
                    help='increase output level')
parser.add_argument('--baseURL', default=baseUrl, help='Base URL for service')
parser.add_argument('--debugLevel', type=int,
                    default=debugLevel, help='debug output level')
options = parser.parse_args()

# Increase output level
if options.verbose:
    outputLevel += 1

# Decrease output level
if options.quiet:
    outputLevel -= 1

# Debug level
if options.debugLevel:
    debugLevel = options.debugLevel

# Debug print


def printDebugMessage(functionName, message, level):
    if(level <= debugLevel):
        print('[' + functionName + '] ' + message, file=sys.stderr)

# User-agent for request (see RFC2616).


def getUserAgent():
    printDebugMessage('getUserAgent', 'Begin', 11)
    # Agent string for urllib2 library.
    urllib_agent = 'Python-urllib/%s' % urllib2.__version__
    clientRevision = '$Revision: 2809 $'
    clientVersion = '0'
    if len(clientRevision) > 11:
        clientVersion = clientRevision[11:-2]
    # Prepend client specific agent string.
    user_agent = 'EBI-Sample-Client/%s (%s; Python %s; %s) %s' % (
        clientVersion, os.path.basename(__file__),
        platform.python_version(), platform.system(),
        urllib_agent
    )
    printDebugMessage('getUserAgent', 'user_agent: ' + user_agent, 12)
    printDebugMessage('getUserAgent', 'End', 11)
    return user_agent

# Wrapper for a REST (HTTP GET) request


def restRequest(url):
    printDebugMessage('restRequest', 'Begin', 11)
    printDebugMessage('restRequest', 'url: ' + url, 11)
    # Errors are indicated by HTTP status codes.
    try:
        # Set the User-agent.
        user_agent = getUserAgent()
        http_headers = {'User-Agent': user_agent}
        req = urllib.request.Request(url, None, http_headers)
        # Make the request (HTTP GET).
        reqH = urllib.request.urlopen(req)
        result = reqH.read()
        reqH.close()
    # Errors are indicated by HTTP status codes.
    except urllib.error.HTTPError as ex:
        # Trap exception and output the document to get error message.
        print(ex.read(), file=sys.stderr)
        raise
    printDebugMessage('restRequest', 'End', 11)
    return result

# Get input parameters list


def serviceGetParameters():
    printDebugMessage('serviceGetParameters', 'Begin', 1)
    requestUrl = baseUrl + '/parameters'
    printDebugMessage('serviceGetParameters', 'requestUrl: ' + requestUrl, 2)
    xmlDoc = restRequest(requestUrl)
    doc = xmltramp.parse(xmlDoc)
    printDebugMessage('serviceGetParameters', 'End', 1)
    return doc['id':]

# Print list of parameters


def printGetParameters():
    printDebugMessage('printGetParameters', 'Begin', 1)
    idList = serviceGetParameters()
    for id in idList:
        print(id)
    printDebugMessage('printGetParameters', 'End', 1)

# Get input parameter information


def serviceGetParameterDetails(paramName):
    printDebugMessage('serviceGetParameterDetails', 'Begin', 1)
    printDebugMessage('serviceGetParameterDetails',
                      'paramName: ' + paramName, 2)
    requestUrl = baseUrl + '/parameterdetails/' + paramName
    printDebugMessage('serviceGetParameterDetails',
                      'requestUrl: ' + requestUrl, 2)
    xmlDoc = restRequest(requestUrl)
    doc = xmltramp.parse(xmlDoc)
    printDebugMessage('serviceGetParameterDetails', 'End', 1)
    return doc

# Print description of a parameter


def printGetParameterDetails(paramName):
    printDebugMessage('printGetParameterDetails', 'Begin', 1)
    doc = serviceGetParameterDetails(paramName)
    print(str(doc.name) + "\t" + str(doc.type))
    print(doc.description)
    for value in doc.values:
        print(value.value, end=' ')
        if str(value.defaultValue) == 'true':
            print('default', end=' ')
        print()
        print("\t" + str(value.label))
        if(hasattr(value, 'properties')):
            for wsProperty in value.properties:
                print("\t" + str(wsProperty.key) + "\t" + str(wsProperty.value))
    #print doc
    printDebugMessage('printGetParameterDetails', 'End', 1)

# Submit job


def serviceRun(email, title, params):
    printDebugMessage('serviceRun', 'Begin', 1)
    # Insert e-mail and title into params
    params['email'] = email
    if title:
        params['title'] = title
    requestUrl = baseUrl + '/run/'
    printDebugMessage('serviceRun', 'requestUrl: ' + requestUrl, 2)
    # Signature methods requires special handling (list)
    applData = ''
    if 'appl' in params:
        # So extract from params
        applList = params['appl']
        del params['appl']
        # Build the method data options
        for appl in applList:
            applData += '&appl=' + appl
    # Get the data for the other options
    requestData = urllib.parse.urlencode(params)
    # Concatenate the two parts.
    requestData += applData
    printDebugMessage('serviceRun', 'requestData: ' + requestData, 2)
    # Errors are indicated by HTTP status codes.
    try:
        # Set the HTTP User-agent.
        user_agent = getUserAgent()
        http_headers = {'User-Agent': user_agent}
        req = urllib.request.Request(requestUrl, None, http_headers)
        # Make the submission (HTTP POST).
        reqH = urllib.request.urlopen(req, requestData)
        jobId = reqH.read()
        reqH.close()
    except urllib.error.HTTPError as ex:
        # Trap exception and output the document to get error message.
        print(ex.read(), file=sys.stderr)
        raise
    printDebugMessage('serviceRun', 'jobId: ' + jobId, 2)
    printDebugMessage('serviceRun', 'End', 1)
    return jobId

# Get job status


def serviceGetStatus(jobId):
    printDebugMessage('serviceGetStatus', 'Begin', 1)
    printDebugMessage('serviceGetStatus', 'jobId: ' + jobId, 2)
    requestUrl = baseUrl + '/status/' + jobId
    printDebugMessage('serviceGetStatus', 'requestUrl: ' + requestUrl, 2)
    status = restRequest(requestUrl)
    printDebugMessage('serviceGetStatus', 'status: ' + status, 2)
    printDebugMessage('serviceGetStatus', 'End', 1)
    return status

# Print the status of a job


def printGetStatus(jobId):
    printDebugMessage('printGetStatus', 'Begin', 1)
    status = serviceGetStatus(jobId)
    print(status)
    printDebugMessage('printGetStatus', 'End', 1)


# Get available result types for job
def serviceGetResultTypes(jobId):
    printDebugMessage('serviceGetResultTypes', 'Begin', 1)
    printDebugMessage('serviceGetResultTypes', 'jobId: ' + jobId, 2)
    requestUrl = baseUrl + '/resulttypes/' + jobId
    printDebugMessage('serviceGetResultTypes', 'requestUrl: ' + requestUrl, 2)
    xmlDoc = restRequest(requestUrl)
    doc = xmltramp.parse(xmlDoc)
    printDebugMessage('serviceGetResultTypes', 'End', 1)
    return doc['type':]

# Print list of available result types for a job.


def printGetResultTypes(jobId):
    printDebugMessage('printGetResultTypes', 'Begin', 1)
    resultTypeList = serviceGetResultTypes(jobId)
    for resultType in resultTypeList:
        print(resultType['identifier'])
        if(hasattr(resultType, 'label')):
            print("\t", resultType['label'])
        if(hasattr(resultType, 'description')):
            print("\t", resultType['description'])
        if(hasattr(resultType, 'mediaType')):
            print("\t", resultType['mediaType'])
        if(hasattr(resultType, 'fileSuffix')):
            print("\t", resultType['fileSuffix'])
    printDebugMessage('printGetResultTypes', 'End', 1)

# Get result


def serviceGetResult(jobId, type_):
    printDebugMessage('serviceGetResult', 'Begin', 1)
    printDebugMessage('serviceGetResult', 'jobId: ' + jobId, 2)
    printDebugMessage('serviceGetResult', 'type_: ' + type_, 2)
    requestUrl = baseUrl + '/result/' + jobId + '/' + type_
    result = restRequest(requestUrl)
    printDebugMessage('serviceGetResult', 'End', 1)
    return result

# Client-side poll


def clientPoll(jobId):
    printDebugMessage('clientPoll', 'Begin', 1)
    result = 'PENDING'
    while result == 'RUNNING' or result == 'PENDING':
        result = serviceGetStatus(jobId)
        print(result, file=sys.stderr)
        if result == 'RUNNING' or result == 'PENDING':
            time.sleep(checkInterval)
    printDebugMessage('clientPoll', 'End', 1)

# Get result for a jobid


def getResult(jobId):
    printDebugMessage('getResult', 'Begin', 1)
    printDebugMessage('getResult', 'jobId: ' + jobId, 1)
    # Check status and wait if necessary
    clientPoll(jobId)
    # Get available result types
    resultTypes = serviceGetResultTypes(jobId)
    for resultType in resultTypes:
        # Derive the filename for the result
        if options.outfile:
            filename = options.outfile + '.' + \
                str(resultType['identifier']) + '.' + \
                str(resultType['fileSuffix'])
        else:
            filename = jobId + '.' + \
                str(resultType['identifier']) + '.' + \
                str(resultType['fileSuffix'])
        # Write a result file
        if not options.outformat or options.outformat == str(resultType['identifier']):
            # Get the result
            result = serviceGetResult(jobId, str(resultType['identifier']))
            fh = open(filename, 'w')
            fh.write(result)
            fh.close()
            print(filename)
    printDebugMessage('getResult', 'End', 1)

# Read a file


def readFile(filename):
    printDebugMessage('readFile', 'Begin', 1)
    fh = open(filename, 'r')
    data = fh.read()
    fh.close()
    printDebugMessage('readFile', 'End', 1)
    return data


# No options... print help.
if numOpts < 2:
    parser.print_help()
# List parameters
elif options.params:
    printGetParameters()
# Get parameter details
elif options.paramDetail:
    printGetParameterDetails(options.paramDetail)
# Submit job
elif options.email and not options.jobid:
    params = {}
    if 1 > 0:
        if os.access(options.input, os.R_OK):  # Read file into content
            params['sequence'] = readFile(options.input)
        else:  # Argument is a sequence id
            params['sequence'] = options.input
    elif options.sequence:  # Specified via option
        if os.access(options.sequence, os.R_OK):  # Read file into content
            params['sequence'] = readFile(options.sequence)
        else:  # Argument is a sequence id
            params['sequence'] = options.sequence
    # Map flag options to boolean values.
    # if options.crc:
    #    params['crc'] = True
    # elif options.nocrc:
    #    params['crc'] = False
    if options.goterms:
        params['goterms'] = True
    elif options.nogoterms:
        params['goterms'] = False
    if options.pathways:
        params['pathways'] = True
    elif options.nopathways:
        params['pathways'] = False
    # Add the other options (if defined)
    if options.appl:
        params['appl'] = re.split('[ \t\n,;]+', options.appl)

    # Submit the job
    jobid = serviceRun(options.email, options.title, params)
    if options.async:  # Async mode
        print(jobid)
    else:  # Sync mode
        print(jobid, file=sys.stderr)
        time.sleep(5)
        getResult(jobid)
# Get job status
elif options.status and options.jobid:
    printGetStatus(options.jobid)
# List result types for job
elif options.resultTypes and options.jobid:
    printGetResultTypes(options.jobid)
# Get results for job
elif options.polljob and options.jobid:
    getResult(options.jobid)
else:
    print('Error: unrecognised argument combination', file=sys.stderr)
    parser.print_help()
