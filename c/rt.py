import resource,subprocess,sys,time
t=time.time(); subprocess.run(sys.argv[1:]); print("RSS_MB=%.0f wall=%.1f"%(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss/1024,time.time()-t),flush=True)
