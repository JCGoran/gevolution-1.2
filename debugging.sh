#!/usr/bin/env bash

# this bash script allows for the use of the GNU debugger (gdb)
# on the root process of gevolution.
# Can be used both locally (when running via `mpirun` or similar),
# or submitted on a remote node (for instance via `qsub` or `sbatch`).
# Note that the remote node is assumed to share the filesystem with
# the submission (frontend) node.
# Upon exit from gdb, gevolution resumes normal operation,
# so if one wishes to terminate the job, one can use `signal SIGTERM`
# while inside the gdb interpreter to terminate the root
# (and correspondingly all the other) process.

file="${HOME}/.__gevolution_debug"
# check the file actually exists
if [ -e "${file}" ]
then
    file_content="$(cat ${file})"
    arr=( ${file_content} )
    p=${arr[0]}
    h=${arr[1]}
    sourcefile=${arr[2]}
    # go to the line in the code stated in the file
    line="$((${arr[3]} + 1))"
    cmd="gdb -p ${p} -ex 'up 2' -ex 'set var i = 1' -ex 'break ${sourcefile}:${line}' -ex 'continue'"
    # if the hostname matches the node the user is currently in, run gdb
    if [ "${h}" = "$(hostname)" ]
    then
        # if PID read from file doesn't exist as a current process, exit
        if [ -z "$(ps h -p ${p})" ]
        then
            printf "ERROR: PID '${p}' doesn't exist\n" 1>&2
            exit 1
        else
            # check if gdb is available
            command -v gdb >/dev/null
            if [ $? -eq 0 ]
            then
                eval "${cmd}"
            else
                printf "ERROR: unable to find 'gdb' on the local system\n" 1>&2
                exit 2
            fi
        fi
    else
        # attempt to ssh to remote
        ssh -q "${h}" exit >/dev/null
        # if it fails for some reason (auth or otherwise), exit
        if [ $? -ne 0 ]
        then
            printf "ERROR: remote '${h}' refused connection\n" 1>&2
            exit 3
        else
            # check if remote has gdb
            ssh -q "${h}" "command -v gdb" >/dev/null
            if [ $? -ne 0 ]
            then
                printf "ERROR: unable to find 'gdb' on remote '${h}'\n" 1>&2
                exit 4
            else
                ssh -t "${h}" "${cmd}"
            fi
        fi
    fi
else
    printf "ERROR: gevolution debug file '${file}' doesn't exist\n" 1>&2
    exit 5
fi
