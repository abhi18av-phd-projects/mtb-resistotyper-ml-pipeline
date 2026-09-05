"""Composable feature-engineering steps.

Every step takes a mart and returns a mart, so techniques compose in any order
the operator asks for without the pipeline knowing what any of them do.

The one thing the pipeline DOES know about a step is whether it reads the
phenotype label. That single property decides where a step may run, and the
registry makes it structural rather than a comment somebody has to remember.
"""
