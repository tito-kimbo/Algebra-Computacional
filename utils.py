class AssumptionnError(Exception):
    pass


"""
    Alternative for asserts
    Usage:
        def myfun(myarg):
            assuming(myarg > 0, "MyArg must be >0")
"""
def assuming(val: bool, msg:str = None):
    if not val:
        if msg is not None:
            raise AssertionError(f"Unfulfilled assumption: {msg}")
        else:
            raise AssertionError(f"Unfulfilled assumption")

