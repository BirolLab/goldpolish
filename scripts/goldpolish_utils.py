import base64
import uuid
import threading
import os
from os.path import join
import signal
import time
import fcntl

INIT_PID = 1
PARENT_QUERY_PERIOD = 1  # In seconds
SIMULTANEOUS_BATCH_PROCESSES_FILE = "simultaneous_batch_processes"
THREADS_IN_USE_FILE = "threads_in_use"


def get_random_name():
    """Return a unique string of characters."""
    return (
        base64.urlsafe_b64encode(uuid.uuid4().bytes)
        .rstrip(b"=")
        .replace(b"-", b"")
        .decode("ascii")
    )


def read_pipe_message(pipepath):
    """Read one newline-terminated message from a named pipe.

    Stops at the newline instead of reading to end of file: on macOS, the final
    end-of-file read() on a named pipe can block forever after the writer has
    closed it. Raw os.read() is used because buffered file reads keep reading
    until end of file.
    """
    fd = os.open(pipepath, os.O_RDONLY)
    try:
        data = b""
        while not data.endswith(b"\n"):
            chunk = os.read(fd, 4096)
            if not chunk:  # Writer closed without a newline.
                break
            data += chunk
    finally:
        os.close(fd)
    return data.decode()


def process_terminate():
    """SIGTERM the calling process."""
    os.kill(os.getpid(), signal.SIGTERM)


def watch_process(process):
    """SIGTERM the calling process if the observed process dies with an error."""

    def _watch_process(process):
        process.wait()
        if process.returncode != 0:
            print(f"{process.args} failed!")
            process_terminate()

    threading.Thread(target=_watch_process, args=(process,), daemon=True).start()


def wait_till_parent_ends():
    """Periodically poll for parent's end."""
    while os.getppid() != INIT_PID:
        time.sleep(PARENT_QUERY_PERIOD)


def bind_to_parent():
    """SIGTERM the calling process if the parent process ends."""

    def _bind_to_parent():
        wait_till_parent_ends()
        process_terminate()

    threading.Thread(target=_bind_to_parent, daemon=True).start()


def create_threadsafe_counter(filepath, initial_value=0):
    lock_path = f"{filepath}.lock"

    with open(filepath, "w") as f:
        print(initial_value, file=f)
    with open(lock_path, "w") as f:
        pass


def get_threadsafe_counter(filepath):
    lock_path = f"{filepath}.lock"
    lock = os.open(lock_path, os.O_RDONLY)

    fcntl.flock(lock, fcntl.LOCK_SH)
    with open(filepath, "r") as counter:
        count = int(counter.read())
    fcntl.flock(lock, fcntl.LOCK_UN)

    os.close(lock)

    return count


def update_threadsafe_counter(filepath, delta):
    lock_path = f"{filepath}.lock"
    lock = os.open(lock_path, os.O_RDONLY)

    fcntl.flock(lock, fcntl.LOCK_EX)
    with open(filepath, "r") as counter:
        count = int(counter.read())
    with open(filepath, "w") as counter:
        print(count + delta, file=counter)
    fcntl.flock(lock, fcntl.LOCK_UN)

    os.close(lock)


def compare_update_threadsafe_counter(filepath, max, delta):
    """
    If the counter value is less than or equal to max - delta,
    then update the counter value and return True.
    """
    lock_path = f"{filepath}.lock"
    lock = os.open(lock_path, os.O_RDONLY)
    updated = False

    fcntl.flock(lock, fcntl.LOCK_EX)
    with open(filepath, "r") as counter:
        count = int(counter.read())
    if count <= max - delta:
        with open(filepath, "w") as counter:
            print(count + delta, file=counter)
        updated = True
    fcntl.flock(lock, fcntl.LOCK_UN)

    os.close(lock)

    return updated


def create_threads_in_use(path, prefix, initial_value):
    threads_in_use_path = join(path, f"{prefix}-{THREADS_IN_USE_FILE}")
    create_threadsafe_counter(threads_in_use_path, initial_value)


def get_threads_in_use(path, prefix):
    threads_in_use_path = join(path, f"{prefix}-{THREADS_IN_USE_FILE}")
    return get_threadsafe_counter(threads_in_use_path)


def update_threads_in_use(path, prefix, delta):
    threads_in_use_path = join(path, f"{prefix}-{THREADS_IN_USE_FILE}")
    update_threadsafe_counter(threads_in_use_path, delta)


def compare_update_threads_in_use(path, prefix, max, delta):
    threads_in_use_path = join(path, f"{prefix}-{THREADS_IN_USE_FILE}")
    return compare_update_threadsafe_counter(threads_in_use_path, max, delta)


def create_simultaneous_batch_processes(path, prefix):
    simultaneous_batch_processes_path = join(
        path, f"{prefix}-{SIMULTANEOUS_BATCH_PROCESSES_FILE}"
    )
    create_threadsafe_counter(simultaneous_batch_processes_path)


def get_simultaneous_batch_processes(path, prefix):
    simultaneous_batch_processes_path = join(
        path, f"{prefix}-{SIMULTANEOUS_BATCH_PROCESSES_FILE}"
    )
    return get_threadsafe_counter(simultaneous_batch_processes_path)


def update_simultaneous_batch_processes(path, prefix, delta):
    simultaneous_batch_processes_path = join(
        path, f"{prefix}-{SIMULTANEOUS_BATCH_PROCESSES_FILE}"
    )
    update_threadsafe_counter(simultaneous_batch_processes_path, delta)
