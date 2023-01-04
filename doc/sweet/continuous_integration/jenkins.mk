# DEPRECATED

This is a remainder document from the try to use Jenkins CI
We now use Gitlab-CI.




# Information about integration in Jenkins

# Pipeline

## Required plugins:

 * Pipeline
 * Github (NOT only the Git plugin)
 * Parallel test executor

## Other infos

Project "Pipeline" required which processes "Jenkinsfile" from repository

Use bash per default:
Manage Jenkins > Configure System > Shell > Shell executable > /bin/bash


# Cache

We want to cache
```
local_software/local
```

https://plugins.jenkins.io/results-cache/

## Required Plugins
 * Job Cacher


# TODO: Parallel build

## Required Plugins

## Further information

https://www.incredibuild.com/blog/jenkins-parallel-builds-jenkins-distributed-builds


# TODO: Docker

## Required Plugins
 * Docker

## Further information

https://www.jenkins.io/doc/book/pipeline/docker





# Checkout repository

git url: 'https://github.com/schreiberx/sweet.git', branch: 'master'

# Fail as soon as one job is faulty
failFast true
