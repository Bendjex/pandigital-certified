"""paths of the C5 code (motor X) and of the certifier package; override with environment variables"""
import os
C5 = os.path.dirname(os.path.abspath(__file__))
SRC = os.environ.get("PANDI_SRC", os.path.join(os.path.dirname(C5), "src"))
RESULTS = os.environ.get("PANDI_RESULTS", os.path.join(os.path.dirname(SRC), "results"))
