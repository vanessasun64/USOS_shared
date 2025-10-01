#Generate directory path for computer type for general saving

import pandas as pd

def filepath_source(computer: str) -> str:
    paths = {
        'CHPC': '/uufs/chpc.utah.edu/common/home/haskins-group1/users/vsun/USOS_shared/',
        'Mac': '/Users/vanessasun/Documents/phd/utah/research/USOS_shared/'
    }

    if computer not in paths:
        raise ValueError(f'Unknown computer identifier. \n Use identifiers "CHPC" or "Mac"')
    
    return paths[computer]