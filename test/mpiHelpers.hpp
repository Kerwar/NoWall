#ifndef MPI_HELPERS
#define MPI_HELPERS

#include <iostream>
#include <mpi.h>

class MPI_Init_H {
    public:
    /* constructor initialises MPI */
    MPI_Init_H( int argc, char* argv[] ){
        MPI_Init(&argc, &argv);
    }
    /* destructor finalises MPI */
    ~MPI_Init_H(){ MPI_Finalize(); }
};

/* content of mpiH_Comm.hpp */
class MPI_Comm_H
{
private:
    MPI_Comm m_comm;
    int m_size;
    int m_rank;

public:
    MPI_Comm_H( MPI_Comm comm = MPI_COMM_NULL ) :
        m_comm {comm}
    {
        if ( m_comm != MPI_COMM_NULL ){
            MPI_Comm_size(m_comm, &m_size);
            MPI_Comm_rank(m_comm, &m_rank);
        }
        else {
            m_size = 0;
            m_rank = -1;
        }
    }

    /* contextual conversion to bool, which returns true if m_comm is a valid
     * communicator and false if it is MPI_COMM_NULL */
    operator bool() const { return m_comm != MPI_COMM_NULL; }
    
    const MPI_Comm& comm() const {
        #ifndef NDEBUG
        if ( !(*this) )
            std::cerr
                << "WARNING: You called comm() on a null communicator!\n";
        #endif
        return m_comm;
    }

    int rank() const {
        return m_rank;
    }

    int size() const {
        assert( *this && "You called size() on a null communicator!");
        return m_size;
    }

    int master() const {
        assert( *this && "You called master() on a null communicator!");
        return 0;
    }

    bool isMaster() const { return rank() == 0; }

    /* allow conversion to MPI_Comm */
    operator const MPI_Comm&() const { return m_comm; }
};

#endif