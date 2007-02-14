#ifndef __NBC_H__
#define __NBC_H__

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <string.h>

#define DLEVEL 3

/* restore inline behavior for non-gcc compilers */
#ifndef __GNUC__
#define __inline__ inline
#endif

/* log(2) */
#define LOG2 0.69314718055994530941

#ifdef __cplusplus
extern "C" {
#endif

/* Function return codes  */
#define NBC_OK 0 /* everything went fine */
#define NBC_SUCCESS 0 /* everything went fine (MPI compliant :) */
#define NBC_OOR 1 /* out of resources */
#define NBC_BAD_SCHED 2 /* bad schedule */
#define NBC_CONTINUE 3 /* progress not done */
#define NBC_DATATYPE_NOT_SUPPORTED 4 /* datatype not supported or not valid */
#define NBC_OP_NOT_SUPPORTED 5 /* operation not supported or not valid */
#define NBC_NOT_IMPLEMENTED 6
#define NBC_INVALID_PARAM 7 /* invalid parameters */

/* dirty trick to avoid fn call :) */
#define NBC_Test NBC_Progress

/* several typedefs for NBC */

/* a schedule is basically a pointer to some memory location where the
 * schedule array resides */ 
typedef void* NBC_Schedule;

/* the function type enum */
typedef enum {
  SEND,
  RECV,
  OP,
  COPY
} NBC_Fn_type;

/* the send argument struct */
typedef struct {
  void *buf;
  int count;
  MPI_Datatype datatype;
  int dest;
} NBC_Args_send;

/* the receive argument struct */
typedef struct {
  void *buf;
  int count;
  MPI_Datatype datatype;
  int source;
} NBC_Args_recv;

/* the operation argument struct */
typedef struct {
  void *buf1;
  void *buf2;
  void *buf3;
  int count;
  MPI_Op op;
  MPI_Datatype datatype;
} NBC_Args_op;

typedef struct {
  void *src; 
  int srccount;
  MPI_Datatype srctype;
  void *tgt;
  int tgtcount;
  MPI_Datatype tgttype;
} NBC_Args_copy;

/* thread specific data */
typedef struct {
  MPI_Comm comm;
  MPI_Comm mycomm;
  long row_offset;
  int tag;
  int req_count;
  MPI_Request *req_array;
  NBC_Schedule *schedule;
  void *tmpbuf; /* temporary buffer e.g. used for Reduce */
} NBC_State;

/* we should make a handle pointer to a state later */
typedef NBC_State NBC_Handle;

/* used to hang off a communicator */
typedef struct {
  MPI_Comm mycomm; /* save the shadow communicator here */
  int tag;
} NBC_Comminfo;

/* internal function prototypes */
int NBC_Sched_create(NBC_Schedule* schedule);
int NBC_Free(NBC_Handle *handle);
int NBC_Sched_send(void* buf, int count, MPI_Datatype datatype, int dest, NBC_Schedule *schedule);
int NBC_Sched_recv(void* buf, int count, MPI_Datatype datatype, int source, NBC_Schedule *schedule);
int NBC_Sched_op(void* buf3, void* buf1, void* buf2, int count, MPI_Datatype datatype, MPI_Op op, NBC_Schedule *schedule);
int NBC_Sched_copy(void *src, int srccount, MPI_Datatype srctype, void *tgt, int tgtcount, MPI_Datatype tgttype, NBC_Schedule *schedule);
int NBC_Sched_barrier(NBC_Schedule *schedule);
int NBC_Sched_commit(NBC_Schedule *schedule);
int NBC_Progress(NBC_Handle *handle);
int NBC_Progress_block(NBC_Handle *handle);
int NBC_Start(NBC_Handle *handle, MPI_Comm comm, NBC_Schedule *schedule);
int NBC_Operation(void *buf3, void *buf1, void *buf2, MPI_Op op, MPI_Datatype type, int count);
static __inline__ int NBC_Type_intrinsic(MPI_Datatype type);
static __inline__ int NBC_Copy(void *src, int srccount, MPI_Datatype srctype, void *tgt, int tgtcount, MPI_Datatype tgttype, MPI_Comm comm);

/* external function prototypes */
int NBC_Ibarrier(MPI_Comm comm, NBC_Handle* handle);
int NBC_Ibcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm, NBC_Handle* handle);
int NBC_Ibcast_lin(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm, NBC_Handle* handle);
int NBC_Ialltoallv(void* sendbuf, int *sendcounts, int *sdispls, MPI_Datatype sendtype, void* recvbuf, int *recvcounts, int *rdispls, MPI_Datatype recvtype, MPI_Comm comm, NBC_Handle* handle);
int NBC_Igather(void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm, NBC_Handle* handle);
int NBC_Igatherv(void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int *recvcounts, int *displs, MPI_Datatype recvtype, int root, MPI_Comm comm, NBC_Handle* handle);
int NBC_Iscatter(void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm, NBC_Handle* handle);
int NBC_Iscatterv(void* sendbuf, int *sendcounts, int *displs, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm, NBC_Handle* handle);
int NBC_Iallgather(void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm, NBC_Handle *handle);
int NBC_Iallgatherv(void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int *recvcounts, int *displs, MPI_Datatype recvtype, MPI_Comm comm, NBC_Handle *handle);
int NBC_Ialltoall(void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm, NBC_Handle *handle);
int NBC_Ireduce(void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm, NBC_Handle* handle);
int NBC_Iallreduce(void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, NBC_Handle* handle);
int NBC_Iscan(void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, NBC_Handle* handle);
int NBC_Ireduce_scatter(void* sendbuf, void* recvbuf, int *recvcounts, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, NBC_Handle* handle);
int NBC_Wait(NBC_Handle *handle);
int NBC_Wait_poll(NBC_Handle *handle);


/* some macros */

/* a schedule has the following format:
 * [schedule] ::= [size][round-schedule][delimiter][round-schedule][delimiter]...[end]
 * [size] ::= size of the schedule (int)
 * [round-schedule] ::= [num][type][type-args][type][type-args]...
 * [num] ::= number of elements in round (int)
 * [type] ::= function type (NBC_Fn_type)
 * [type-args] ::= type specific arguments (NBC_Args_send, NBC_Args_recv or, NBC_Args_op)
 * [delimiter] ::= 1 (char) - indicates that a round follows
 * [end] ::= 0 (char) - indicates that this is the last round 
 */

/* NBC_GET_ROUND_SIZE returns the size in bytes of a round of a NBC_Schedule
 * schedule. A round has the format:
 * [num]{[type][type-args]}
 * e.g. [(int)2][(NBC_Fn_type)SEND][(NBC_Args_send)SEND-ARGS][(NBC_Fn_type)RECV][(NBC_Args_recv)RECV-ARGS] */
#define NBC_GET_ROUND_SIZE(schedule, size) \
 {  \
   int *numptr; \
   NBC_Fn_type *typeptr; \
   int i;  \
     \
   numptr = (int*)schedule; \
   /*DEBUG(10, "GET_ROUND_SIZE got %i elements\n", *numptr); */\
   /* end is increased by sizeof(int) bytes to point to type */ \
   typeptr = (NBC_Fn_type*)((int*)(schedule)+1); \
   for (i=0; i<*numptr; i++) { \
     /* go sizeof op-data forward */ \
     switch(*typeptr) { \
       case SEND: \
         /*printf("found a SEND at offset %i\n", (int)typeptr-(int)schedule); */\
         typeptr = (NBC_Fn_type*)((NBC_Args_send*)typeptr+1); \
         break; \
       case RECV: \
         /*printf("found a RECV at offset %i\n", (int)typeptr-(int)schedule); */\
         typeptr = (NBC_Fn_type*)((NBC_Args_recv*)typeptr+1); \
         break; \
       case OP: \
         /*printf("found a OP at offset %i\n", (int)typeptr-(int)schedule); */\
         typeptr = (NBC_Fn_type*)((NBC_Args_op*)typeptr+1); \
         break; \
       case COPY: \
         /*printf("found a COPY at offset %i\n", (int)typeptr-(int)schedule); */\
         typeptr = (NBC_Fn_type*)((NBC_Args_copy*)typeptr+1); \
         break; \
       default: \
         printf("NBC_GET_ROUND_SIZE: bad type %li at offset %li\n", (long)*typeptr, (long)typeptr-(long)schedule); \
         return NBC_BAD_SCHED; \
     } \
     /* increase ptr by size of fn_type enum */ \
     typeptr = (NBC_Fn_type*)((NBC_Fn_type*)typeptr+1); \
   } \
   /* this could be optimized if typeptr would be used directly */ \
   size = (long)typeptr-(long)schedule; \
 }

/* returns the size of a schedule in bytes */
#define NBC_GET_SIZE(schedule, size) \
{ \
  size=*(int*)schedule; \
}

/* increase the size of a schedule by size bytes */
#define NBC_INC_SIZE(schedule, size) \
{ \
  *(int*)schedule+=size; \
}

/* increments the number of operations in the last round */
#define NBC_INC_NUM_ROUND(schedule) \
{ \
  int total_size; \
  long round_size; \
  char *ptr, *lastround; \
 \
  NBC_GET_SIZE(schedule, total_size); \
 \
  /* ptr begins at first round (first int is overall size) */ \
  ptr = (char*)((char*)schedule+sizeof(int)); \
  lastround = ptr; \
  while ((long)ptr-(long)schedule < total_size) { \
    NBC_GET_ROUND_SIZE(ptr, round_size); \
    /*printf("got round size %i\n", round_size);*/ \
    lastround = ptr; \
    /* add round size */ \
    ptr=ptr+round_size; \
    /* add sizeof(char) as barrier delimiter */ \
    ptr=ptr+sizeof(char); \
    /*printf("(int)ptr-(int)schedule=%i, size=%i\n", (int)ptr-(int)schedule, size); */\
  } \
  /*printf("lastround count is at offset: %i\n", (int)lastround-(int)schedule);*/ \
  /* this is the count in the last round of the schedule */ \
  (*(int*)lastround)++; \
}

/* NBC_PRINT_ROUND prints a round in a schedule. A round has the format:
 * [num]{[op][op-data]} types: [int]{[enum][op-type]}
 * e.g. [2][SEND][SEND-ARGS][RECV][RECV-ARGS] */
#define NBC_PRINT_ROUND(schedule) \
 {  \
   int myrank, *numptr; \
   NBC_Fn_type *typeptr; \
   NBC_Args_send *sendargs; \
   NBC_Args_recv *recvargs; \
   NBC_Args_op *opargs; \
   NBC_Args_copy *copyargs; \
   int i;  \
     \
   numptr = (int*)schedule; \
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank); \
   printf("has %i actions: \n", *numptr); \
   /* end is increased by sizeof(int) bytes to point to type */ \
   typeptr = (NBC_Fn_type*)((int*)(schedule)+1); \
   for (i=0; i<*numptr; i++) { \
     /* go sizeof op-data forward */ \
     switch(*typeptr) { \
       case SEND: \
         printf("[%i]  SEND (offset %li) ", myrank, (long)typeptr-(long)schedule); \
         sendargs = (NBC_Args_send*)(typeptr+1); \
         printf("*buf: %lu, count: %i, type: %lu, dest: %i)\n", (unsigned long)sendargs->buf, sendargs->count, (unsigned long)sendargs->datatype, sendargs->dest); \
         typeptr = (NBC_Fn_type*)((NBC_Args_send*)typeptr+1); \
         break; \
       case RECV: \
         printf("[%i]  RECV (offset %li) ", myrank, (long)typeptr-(long)schedule); \
         recvargs = (NBC_Args_recv*)(typeptr+1); \
         printf("*buf: %lu, count: %i, type: %lu, source: %i)\n", (unsigned long)recvargs->buf, recvargs->count, (unsigned long)recvargs->datatype, recvargs->source); \
         typeptr = (NBC_Fn_type*)((NBC_Args_recv*)typeptr+1); \
         break; \
       case OP: \
         printf("[%i]  OP   (offset %li) ", myrank, (long)typeptr-(long)schedule); \
         opargs = (NBC_Args_op*)(typeptr+1); \
         printf("*buf1: %lu, buf2: %lu, count: %i, type: %lu)\n", (unsigned long)opargs->buf1, (unsigned long)opargs->buf2, opargs->count, (unsigned long)opargs->datatype); \
         typeptr = (NBC_Fn_type*)((NBC_Args_op*)typeptr+1); \
         break; \
       case COPY: \
         printf("[%i]  COPY   (offset %li) ", myrank, (long)typeptr-(long)schedule); \
         copyargs = (NBC_Args_copy*)(typeptr+1); \
         printf("*src: %lu, srccount: %i, srctype: %lu, *tgt: %lu, tgtcount: %i, tgttype: %lu)\n", (unsigned long)copyargs->src, copyargs->srccount, (unsigned long)copyargs->srctype, (unsigned long)copyargs->tgt, copyargs->tgtcount, (unsigned long)copyargs->tgttype); \
         typeptr = (NBC_Fn_type*)((NBC_Args_copy*)typeptr+1); \
         break; \
       default: \
         printf("[%i] NBC_PRINT_ROUND: bad type %li at offset %li\n", myrank, (long)*typeptr, (long)typeptr-(long)schedule); \
         return NBC_BAD_SCHED; \
     } \
     /* increase ptr by size of fn_type enum */ \
     typeptr = (NBC_Fn_type*)((NBC_Fn_type*)typeptr+1); \
   } \
   printf("\n"); \
 }

#define NBC_PRINT_SCHED(schedule) \
{ \
  int size, myrank; \
  long round_size; \
  char *ptr; \
 \
  NBC_GET_SIZE(schedule, size); \
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank); \
  printf("[%i] printing schedule of size %i\n", myrank, size); \
 \
  /* ptr begins at first round (first int is overall size) */ \
  ptr = (char*)((char*)schedule+sizeof(int)); \
  while ((long)ptr-(long)schedule < size) { \
    NBC_GET_ROUND_SIZE(ptr, round_size); \
    printf("[%i] Round at byte %li (size %li) ", myrank, (long)ptr-(long)schedule, round_size); \
    NBC_PRINT_ROUND(ptr); \
    /* add round size */ \
    ptr=ptr+round_size; \
    /* add sizeof(char) as barrier delimiter */ \
    ptr=ptr+sizeof(char); \
  } \
}

#define CHECK_NULL(ptr) \
{ \
  if(ptr == NULL) { \
    printf("realloc error :-(\n"); \
  } \
}

#define DEBUG(level, ...) {} 
/*
static __inline__ void DEBUG(int level, const char *fmt, ...) 
{ 
  va_list ap;
  int rank; 
 
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
   
  if(DLEVEL > level) { 
    printf("[%i] ", rank); 
    printf(fmt, ap);
    va_end (ap);
  } 
}*/

/* returns true (1) or false (0) if type is intrinsic or not */
static __inline__ int NBC_Type_intrinsic(MPI_Datatype type) {
  
  if( ( type == MPI_INT ) ||
      ( type == MPI_LONG ) ||
      ( type == MPI_SHORT ) ||
      ( type == MPI_UNSIGNED ) ||
      ( type == MPI_UNSIGNED_SHORT ) ||
      ( type == MPI_UNSIGNED_LONG ) ||
      ( type == MPI_FLOAT ) ||
      ( type == MPI_DOUBLE ) ||
      ( type == MPI_LONG_DOUBLE ) ||
      ( type == MPI_BYTE ) ||
      ( type == MPI_FLOAT_INT) ||
      ( type == MPI_DOUBLE_INT) ||
      ( type == MPI_LONG_INT) ||
      ( type == MPI_2INT) ||
      ( type == MPI_SHORT_INT) ||
      ( type == MPI_LONG_DOUBLE_INT)) 
    return 1;
  else 
    return 0;
}

/* let's give a try to inline functions */
static __inline__ int NBC_Copy(void *src, int srccount, MPI_Datatype srctype, void *tgt, int tgtcount, MPI_Datatype tgttype, MPI_Comm comm) {
  int size, pos, res;
  MPI_Aint ext;
  void *packbuf;

  if((srctype == tgttype) && NBC_Type_intrinsic(srctype)) {
    /* if we have the same types and they are contiguous (intrinsic
     * types are contiguous), we can just use a single memcpy */
    res = MPI_Type_extent(srctype, &ext);
    if (MPI_SUCCESS != res) { printf("MPI Error in MPI_Type_extent() (%i)\n", res); return res; }
    memcpy(tgt, src, srccount*ext);
  } else {
    /* we have to pack and unpack */
    res = MPI_Pack_size(srccount, srctype, comm, &size);
    if (MPI_SUCCESS != res) { printf("MPI Error in MPI_Pack_size() (%i)\n", res); return res; }
    packbuf = malloc(size);
    if (NULL == packbuf) { printf("Error in malloc()\n"); return res; }
    pos=0;
    res = MPI_Pack(src, srccount, srctype, packbuf, size, &pos, comm);
    if (MPI_SUCCESS != res) { printf("MPI Error in MPI_Pack() (%i)\n", res); return res; }
    pos=0;
    res = MPI_Unpack(packbuf, size, &pos, tgt, tgtcount, tgttype, comm);
    if (MPI_SUCCESS != res) { printf("MPI Error in MPI_Unpack() (%i)\n", res); return res; }
    free(packbuf);
  }

  return NBC_OK;
}

#ifdef __cplusplus
}
#endif
#endif
