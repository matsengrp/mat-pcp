# Contributing to DASM2

This guide covers code style, the PR process, and common pitfalls.

## Development Setup

1. Clone the repository
2. Create virtual environment: `python -m venv .venv`
3. Activate: `source .venv/bin/activate`
4. Install dev dependencies: `pip install -e ".[dev]"`

## Code Style

### Type Hints

All functions and classes must have complete type annotations:

```python
def compute_loss(
    predictions: torch.Tensor,
    targets: torch.Tensor,
    mask: torch.Tensor,
) -> torch.Tensor:
    ...
```

### Docstrings

Use Google-style docstrings:

```python
def optimize_branch_lengths(
    dataset: EvolutionaryDataset,
    model: CoreModel,
) -> Dict[str, float]:
    """Optimize branch lengths for all sequences in dataset.

    Args:
        dataset: Dataset containing sequences and current branch lengths.
        model: Trained model for computing likelihoods.

    Returns:
        Dictionary mapping sequence IDs to optimized branch lengths.

    Raises:
        ValueError: If dataset contains no sequences.
    """
```

### Naming Conventions

- **Count variables**: `num_layers`, `num_sites`, `num_sequences`
- **Boolean variables as questions**: `has_structural_data`, `is_paired`
- **Intent-revealing names**: `evolutionary_loss_weight` not `w1`
- **Collection suffixes**: `csps_list` for lists, `csps_batch` for batched tensors

### Import Organization

```python
# Standard library
import os
from pathlib import Path

# Third-party
import torch
import pandas as pd

# Local
from dasm2.losses import compute_dasm_loss
from dasm2.data import EvolutionaryDataset
```

## Testing

### Real Tests Only

Tests must validate real behavior:

```python
# GOOD: Test against real netam implementation
def test_whichmut_loss_matches_netam():
    result = compute_whichmut_loss(...)
    expected = netam_whichmut_loss(...)
    torch.testing.assert_close(result, expected)

# BAD: Fake mock that doesn't validate anything
def test_model():
    mock_model = Mock(return_value=torch.zeros(10, 20))  # Don't do this
```

### Run Tests

```bash
make test        # Run all tests
make test-cov    # With coverage report
```

## Pre-PR Checklist

Before submitting a pull request:

1. **Format code**: `make format`
2. **Run linter**: `make lint`
3. **Run type checker**: `make mypy`
4. **Run tests**: `make test`
5. **All checks**: `make check`

## Common Pitfalls

### Silent Defaults

Don't use `.get()` with defaults for required configuration:

```python
# BAD: Silent fallback hides configuration problems
d_model = config.get("d_model", 256)

# GOOD: Fail fast with clear error
d_model = config["d_model"]  # KeyError if missing
```

### God Classes

Keep classes focused on one responsibility:

```python
# BAD: Trainer doing everything
class MonolithicTrainer:
    def train(self): ...
    def optimize_branch_lengths(self): ...
    def save_checkpoint(self): ...
    def log_metrics(self): ...

# GOOD: Separate concerns
class Trainer:
    def __init__(self, hooks): ...
    def train(self): ...  # Delegates to hooks

class BranchLengthHook: ...
class CheckpointHook: ...
```

### Hardcoded Switches

Prefer composition over configuration flags:

```python
# BAD: String switches
if head_type == "dasm":
    return DASMHead()
elif head_type == "whichmut":
    return WhichmutHead()

# GOOD: Registry pattern
HEAD_REGISTRY = {"dasm": DASMHead, "whichmut": WhichmutHead}
return HEAD_REGISTRY[head_type]()
```

### Silent Error Handling

Never swallow exceptions:

```python
# BAD: Silent failure
try:
    trainer = load_checkpoint(path)
except:
    trainer = None

# GOOD: Explicit error
try:
    trainer = load_checkpoint(path)
except FileNotFoundError:
    raise ValueError(f"Checkpoint not found at {path}")
```

## PR Guidelines

### Description Format

```markdown
## Summary
Brief description of what changed.

## Changes
- First change
- Second change

Closes #123
```

### Review Process

1. Ensure all CI checks pass
2. Request review from maintainers
3. Address feedback
4. Squash and merge when approved

## Documentation Updates

If your changes affect user-facing behavior, update the relevant docs:

- New features → Update guides/
- API changes → Update reference/
- New concepts → Update concepts/
