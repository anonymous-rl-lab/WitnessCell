"""Public exception hierarchy for WitnessCell."""


class WitnessCellError(Exception):
    """Base class for package-specific failures."""


class ValidationError(WitnessCellError, ValueError):
    """Raised when data or split contracts are invalid."""


class NotFittedError(WitnessCellError, RuntimeError):
    """Raised when prediction is requested before fitting."""


class SerializationError(WitnessCellError, RuntimeError):
    """Raised for unsafe, corrupted, or incompatible model bundles."""

