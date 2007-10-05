/*
 * Copyright (c) 2006 The Trustees of Indiana University and Indiana
 *                    University Research and Technology
 *                    Corporation.  All rights reserved.
 * Copyright (c) 2006 The Technical University of Chemnitz. All 
 *                    rights reserved.
 */
#ifndef MCA_COLL_libnbc_EXPORT_H
#define MCA_COLL_libnbc_EXPORT_H

#include "ompi_config.h"

#include "mpi.h"
#include "opal/mca/mca.h"
#include "ompi/mca/coll/coll.h"
#include "ompi/request/request.h"
#include "ompi/mca/pml/pml.h"

#include "nbc.h"

#if defined(c_plusplus) || defined(__cplusplus)
extern "C" {
#endif

/*
 * Globally exported variable
 */

OMPI_COMP_EXPORT extern const mca_coll_base_component_1_0_0_t mca_coll_libnbc_component;
OMPI_COMP_EXPORT extern int mca_coll_libnbc_priority_param;


/*
 * coll API functions
 */


  /* API functions */

  int mca_coll_libnbc_init_query(bool enable_progress_threads,
                               bool enable_mpi_threads);
  const struct mca_coll_base_module_1_0_0_t *
    mca_coll_libnbc_comm_query(struct ompi_communicator_t *comm, int *priority,
                             struct mca_coll_base_comm_t **data);

  const struct mca_coll_base_module_1_0_0_t *
    mca_coll_libnbc_module_init(struct ompi_communicator_t *comm);
  int mca_coll_libnbc_module_finalize(struct ompi_communicator_t *comm);

  int mca_coll_libnbc_allgather_intra(void *sbuf, int scount, 
                                     struct ompi_datatype_t *sdtype, 
                                     void *rbuf, int rcount, 
                                     struct ompi_datatype_t *rdtype, 
                                     struct ompi_communicator_t *comm);
  int mca_coll_libnbc_allgatherv_intra(void *sbuf, int scount, 
                                      struct ompi_datatype_t *sdtype, 
                                      void * rbuf, int *rcounts, int *disps, 
                                      struct ompi_datatype_t *rdtype, 
                                      struct ompi_communicator_t *comm);
  int mca_coll_libnbc_allreduce_intra(void *sbuf, void *rbuf, int count, 
                                     struct ompi_datatype_t *dtype, 
                                     struct ompi_op_t *op, 
                                     struct ompi_communicator_t *comm);
  int mca_coll_libnbc_alltoall_intra(void *sbuf, int scount, 
                                    struct ompi_datatype_t *sdtype, 
                                    void* rbuf, int rcount, 
                                    struct ompi_datatype_t *rdtype, 
                                    struct ompi_communicator_t *comm);
  int mca_coll_libnbc_alltoallv_intra(void *sbuf, int *scounts, int *sdisps, 
                                     struct ompi_datatype_t *sdtype, 
                                     void *rbuf, int *rcounts, int *rdisps, 
                                     struct ompi_datatype_t *rdtype, 
                                     struct ompi_communicator_t *comm);
  int mca_coll_libnbc_alltoallw_intra(void *sbuf, int *scounts, int *sdisps, 
                                     struct ompi_datatype_t **sdtypes, 
                                     void *rbuf, int *rcounts, int *rdisps, 
                                     struct ompi_datatype_t **rdtypes, 
                                     struct ompi_communicator_t *comm);
  int mca_coll_libnbc_barrier_intra(struct ompi_communicator_t *comm);
  int mca_coll_libnbc_bcast_intra(void *buff, int count, 
                                struct ompi_datatype_t *datatype,
                                int root, 
                                struct ompi_communicator_t *comm);
  int mca_coll_libnbc_exscan_intra(void *sbuf, void *rbuf, int count, 
                                  struct ompi_datatype_t *dtype, 
                                  struct ompi_op_t *op, 
                                  struct ompi_communicator_t *comm);
  int mca_coll_libnbc_gather_intra(void *sbuf, int scount, 
                                  struct ompi_datatype_t *sdtype, void *rbuf, 
                                  int rcount, struct ompi_datatype_t *rdtype, 
                                  int root, struct ompi_communicator_t *comm);
  int mca_coll_libnbc_gatherv_intra(void *sbuf, int scount, 
                                   struct ompi_datatype_t *sdtype, void *rbuf, 
                                   int *rcounts, int *disps, 
                                   struct ompi_datatype_t *rdtype, int root, 
                                   struct ompi_communicator_t *comm);
  int mca_coll_libnbc_reduce_intra(void *sbuf, void* rbuf, int count, 
                                 struct ompi_datatype_t *dtype, 
                                 struct ompi_op_t *op, 
                                 int root,
                                 struct ompi_communicator_t *comm);
  int mca_coll_libnbc_reduce_scatter_intra(void *sbuf, void *rbuf, 
                                          int *rcounts, 
                                          struct ompi_datatype_t *dtype, 
                                          struct ompi_op_t *op, 
                                          struct ompi_communicator_t *comm);
  int mca_coll_libnbc_scan_intra(void *sbuf, void *rbuf, int count, 
                                struct ompi_datatype_t *dtype, 
                                struct ompi_op_t *op, 
                                struct ompi_communicator_t *comm);
  int mca_coll_libnbc_scatter_intra(void *sbuf, int scount, 
                                   struct ompi_datatype_t *sdtype, void *rbuf, 
                                   int rcount, struct ompi_datatype_t *rdtype, 
                                   int root, struct ompi_communicator_t *comm);
  int mca_coll_libnbc_scatterv_intra(void *sbuf, int *scounts, int *disps, 
                                    struct ompi_datatype_t *sdtype, 
                                    void* rbuf, int rcount, 
                                    struct ompi_datatype_t *rdtype, int root, 
                                    struct ompi_communicator_t *comm);

#if defined(c_plusplus) || defined(__cplusplus)
}
#endif
#endif /* MCA_COLL_libnbc_EXPORT_H */
