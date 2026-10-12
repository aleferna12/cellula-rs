# Declaration

The code in this repository was entirely conceptualised by a human, without the assistance of AI agents and tools.

AI was used exclusively for:

1. Writing documentation
2. Writing tests
3. Formatting code
4. Single-line code completion

# Motivation

Code meant for use in research must be trustworthy, auditable, and clear. 
Heavy use of AI agents in research software produces code that meets none of these criteria: LLMs are [overconfident in their intellectual capabilities](https://doi.org/10.48550/arXiv.2505.02151) and much less hesitant to make mistakes than humans when facing uncertainty.
As a result, subtle bugs that would normally be caught during human-driven development can easily infiltrate scientific databases and persist unchecked.
Furthermore, [AI writes bad code](https://www.scrum.org/resources/blog/can-ai-generate-high-quality-code), creating [maintainability and clarity challenges](https://unicoconnect.com/blogs/ai-code-at-scale) that stem from the way LLMs are trained to achieve immediate goals while overlooking the more fundamental questions of code design and architecture.
