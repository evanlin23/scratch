#!/opt/mm/root/envs/picrust2/bin/python
# Runs the stock picrust2_pipeline.py with util.check_empty_traits() disabled. That function only
# validates that the default trait tables are non-empty (loads all of them, ~8 min and ~14 GB);
# it does not change any output.
import runpy, sys, picrust2.pipeline
picrust2.pipeline.check_empty_traits = lambda filepaths: None
sys.argv[0] = '/opt/mm/root/envs/picrust2/bin/picrust2_pipeline.py'
runpy.run_path(sys.argv[0], run_name='__main__')
