/* Copyright (c) 2018 Jakob Meng, <jakobmeng@web.de>
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include <hbrs/mpl/detail/mpi.hpp>
#include <hbrs/mpl/core/preprocessor.hpp>
#include <hbrs/mpl/dt/exception.hpp>
#include <boost/throw_exception.hpp>
#include <boost/assert.hpp>
#include <mpi.h>

HBRS_MPL_NAMESPACE_BEGIN
namespace detail {
namespace mpi {

static void
safe(int ec) {
	if (ec != MPI_SUCCESS) {
		BOOST_THROW_EXCEPTION((
			mpi_exception{} 
			<< errinfo_mpi_error_info{mpi_error_info{ec}}
		));
	}
}

MPI_Datatype
datatype(hana::basic_type<double>) {
	return MPI_DOUBLE;
}
MPI_Datatype
datatype(hana::basic_type<int>) {
	return MPI_INT;
}
MPI_Datatype
datatype(hana::basic_type<unsigned long>) {
	return MPI_UNSIGNED_LONG;
}
MPI_Datatype
datatype(hana::basic_type<unsigned long long>) {
	return MPI_UNSIGNED_LONG_LONG;
}

bool
initialized() {
	int flag;
	safe(MPI_Initialized(&flag));
	return flag;
}

int
rank(MPI_Comm comm) {
	BOOST_ASSERT(initialized());
	int rank;
	safe(MPI_Comm_rank(comm, &rank));
	return rank;
}

int
size(MPI_Comm comm) {
	BOOST_ASSERT(initialized());
	int size;
	safe(MPI_Comm_size(comm, &size));
	return size;
}

void
barrier(MPI_Comm comm) {
	BOOST_ASSERT(initialized());
	safe(MPI_Barrier(comm));
}

void
abort(MPI_Comm comm, int errorcode) {
	BOOST_ASSERT(initialized());
	safe(MPI_Abort(comm, errorcode));
}

MPI_Request
ibcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm) {
	BOOST_ASSERT(initialized());
	MPI_Request request;
	safe(MPI_Ibcast(buffer, count, datatype, root, comm, &request));
	return request;
}

MPI_Status
wait(MPI_Request & request) {
	BOOST_ASSERT(initialized());
	MPI_Status status;
	safe(MPI_Wait(&request, &status));
	return status;
}

MPI_Status
probe(int source, int tag, MPI_Comm comm) {
	BOOST_ASSERT(initialized());
	MPI_Status status;
	safe(MPI_Probe(source, tag, comm, &status));
	return status;
}

int
get_count(MPI_Status const& status, MPI_Datatype datatype) {
	BOOST_ASSERT(initialized());
	int count;
	safe(MPI_Get_count(&status, datatype, &count));
	return count;
}

MPI_Request
isend(void const *buffer, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm) {
	BOOST_ASSERT(initialized());
	MPI_Request request;
	safe(MPI_Isend(buffer, count, datatype, dest, tag, comm, &request));
	return request;
}

MPI_Request
irecv(void *buffer, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm) {
	BOOST_ASSERT(initialized());
	MPI_Request request;
	safe(MPI_Irecv(buffer, count, datatype, source, tag, comm, &request));
	return request;
}

/* namespace mpi */ }
/* namespace detail */ }
HBRS_MPL_NAMESPACE_END
