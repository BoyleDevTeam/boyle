# Python Conventions

## Style

- **PEP 8** baseline; **ruff** enforces project-specific rules via `pyproject.toml`.
- **Type hints** on all public functions and methods.
- Imports grouped: stdlib → third party → local.

## Tools

| Tool | Role |
| ---- | ---- |
| `uv` | Environment and lockfile |
| `ruff` | Lint + format |
| `pyright` | Static typing |
| `pytest` | Tests (`uv run pytest`) |

## Modules

- Avoid empty `__init__.py` unless packaging requires it; prefer explicit imports.

## Tests

- **Arrange / Act / Assert** structure
- **`@pytest.fixture`** for reusable setup
- **`@pytest.mark.parametrize`** for data-driven cases

## Adding dependencies

```bash
uv add package_name
uv lock
```

See `.github/instructions/python-coding-standards.instructions.md` for more detail.
