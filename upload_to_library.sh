image=$1
librarypath=$2 #something like usrid/collection/imagename:version

#Eventually login into remote
#singularity remote login --tokenfile singularity_library.token
singularity push -U $image library://$librarypath
