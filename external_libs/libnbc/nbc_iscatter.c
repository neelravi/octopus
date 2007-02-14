/*
 * Copyright (c) 2006 The Trustees of Indiana University and Indiana
 *                    University Research and Technology
 *                    Corporation.  All rights reserved.
 * Copyright (c) 2006 The Technical University of Chemnitz. All 
 *                    rights reserved.
 */
#include "nbc.h"

/* simple linear MPI_Iscatter */
int NBC_Iscatter(void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm, NBC_Handle* handle) {
  int rank, p, res, i;
  MPI_Aint sndext;
  NBC_Schedule *schedule;
  char *sbuf;
  
  res = MPI_Comm_rank(comm, &rank);
  if (MPI_SUCCESS != res) { printf("MPI Error in MPI_Comm_rank() (%i)\n", res); return res; }
  res = MPI_Comm_size(comm, &p);
  if (MPI_SUCCESS != res) { printf("MPI Error in MPI_Comm_size() (%i)\n", res); return res; }
  res = MPI_Type_extent(sendtype, &sndext);
  if (MPI_SUCCESS != res) { printf("MPI Error in MPI_Type_extent() (%i)\n", res); return res; }

  schedule = malloc(sizeof(NBC_Schedule));
  if (NULL == schedule) { printf("Error in malloc()\n"); return res; }

  handle->tmpbuf=NULL;
 
  res = NBC_Sched_create(schedule);
  if(res != NBC_OK) { printf("Error in NBC_Sched_create (%i)\n", res); return res; }

  /* receive from root */
  if(rank != root) {
    /* recv msg from root */
    res = NBC_Sched_recv(recvbuf, recvcount, recvtype, root, schedule);
    if (NBC_OK != res) { printf("Error in NBC_Sched_recv() (%i)\n", res); return res; }
  } else {
    for(i=0;i<p;i++) {
      sbuf = ((char *)sendbuf) + (i*sendcount*sndext);
      if(i == root) {
        /* if I am the root - just copy the message */
        res = NBC_Copy(sbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);
        if (NBC_OK != res) { printf("Error in NBC_Copy() (%i)\n", res); return res; }
      } else {
        /* root sends the right buffer to the right receiver */
        res = NBC_Sched_send(sbuf, sendcount, sendtype, i, schedule);
        if (NBC_OK != res) { printf("Error in NBC_Sched_send() (%i)\n", res); return res; }
      }
    }
  }
 
  res = NBC_Sched_commit(schedule);
  if (NBC_OK != res) { printf("Error in NBC_Sched_commit() (%i)\n", res); return res; }
 
  res = NBC_Start(handle, comm, schedule);
  if (NBC_OK != res) { printf("Error in NBC_Start() (%i)\n", res); return res; }
 
  return NBC_OK;
}
