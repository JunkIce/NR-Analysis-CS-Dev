print('Hello')
import threading
import time
def spin(stop_event):
    i=0
    strs=['|','/','-','\\']
    while not stop_event.is_set():
        print(f'Loading {strs[i%len(strs)]}',end='\r')
        time.sleep(0.15)
        i+=1
    print('         ',end='\r')
    
stop_event=threading.Event()
t1=threading.Thread(target=spin, args=(stop_event,))

t1.start()

time.sleep(5)

stop_event.set()
t1.join()
print('JOINED')