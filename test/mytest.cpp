#include <memory>
#include "dd/Package.hpp"

auto dd = std::make_unique<dd::Package>(1); // Create new package instance capable of handling a single qubit
auto zero_state = dd->makeZeroState(1) ; // zero_state = |0>

/* Creating a DD requires the following inputs:
 * 1. A 2x2 matrix describing a single-qubit operation (here: the Hadamard matrix)
 * 2. The number of qubits the DD will operate on (here: one qubit)
 * 3. The qubit the operation is applied to (here: q0) 
 * (4. Controlled operations can be created by additionally specifying a list of control qubits before the target declaration)
 */
auto h_op = dd->makeGateDD(dd::Hmat, 1, 0);

// Multiplying the operation and the state results in a new state, here a single qubit in superposition
auto superposition = dd->multiply(h_op, zero_state); 

int main()
{
    return 0;
}