## Core Testing Philosophy

Unit tests should validate correctness through systematic, progressive
verification that builds confidence from fundamental operations to complex
behaviors. Every test must serve the purpose of validating the intended
functionality of the target being tested.

## 0. Purpose-Driven Testing

- **Validate intended functionality**: Every test must directly validate a
  specific aspect of the target's designed behavior
- **No superfluous tests**: If a test doesn't verify actual functionality,
  remove it
- **Analyze necessity# Hierarchical and Incremental, Specific and Detailed Unit
  Testing Principles

## Core Testing Philosophy

Unit tests should validate correctness through systematic, progressive
verification that builds confidence from fundamental operations to complex
behaviors. Every test must serve the purpose of validating the intended
functionality of the target being tested.

## 1. Hierarchical Testing

- **Start with fundamentals**: Test basic operations first (object creation,
  simple getters/setters, single operations)
- **Build systematically**: Complex scenarios should only be tested after their
  component parts are validated
- **Isolate failures quickly**: When a complex test fails, simpler tests should
  pinpoint the root cause
- **Example**: For a memory allocator, test single allocation → multiple
  allocations → allocations with alignment → allocations across partitions

## 2. Incremental Validation

- **Progressive complexity**: Add one dimension of complexity at a time
- **No gaps**: Don't jump from trivial to complex cases - test intermediate
  scenarios
- **Combinatorial growth**: When testing multiple features, vary one at a time
  before combining
- **Example**: For distributed arrays, test 1 partition → 2 partitions →
  3,4,5... partitions → variable-sized partitions

## 3. Specific Property Testing

- **Testable assertions**: Every test should verify properties that have
  definite correct values
- **Invariants**: Identify and test properties that must always hold (e.g.,
  allocated size = requested size)
- **Arithmetic correctness**: When numbers are involved, the expected values
  must be arithmetically derivable
- **No arbitrary checks**: Don't test for values without a concrete reason why
  that value is correct
- **Example**: Total memory used = sum of partition sizes; offsets[i+1] -
  offsets[i] = partition[i] size

## 4. Detailed Coverage

- **Edge cases**: Empty inputs, single elements, boundary values
- **Error conditions**: Invalid inputs, resource exhaustion
- **Special configurations**: Test non-default settings and their interactions
- **Real-world scenarios**: Include cases that mirror actual usage patterns

## 5. Implementation Guidelines

- **Test independence**: Each test should set up its own state and not depend on
  test order
- **Clear test names**: Name should indicate what is being tested and expected
  outcome
- **Appropriate tolerances**: Use exact comparisons for discrete values, small
  tolerances only for floating point
- **Progressive test structure**: Order tests within a suite from basic to
  complex
- **Focused assertions**: Each test should validate one specific behavior or
  property

## 6. Error Diagnosis

- **Clear failure messages**: When a test fails, the output should indicate what
  was expected vs actual
- **Intermediate state validation**: Complex operations should check
  intermediate results
- **Minimal reproducible cases**: Complex test failures should be reducible to
  simpler failing tests

## Key Principle

**"Test the simplest thing that could possibly fail"** - Start with basic
functionality and systematically add complexity. Each test should have a clear
reason for existing based on a specific property or behavior that must be
verified.