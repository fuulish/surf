import os
import pathlib
import shutil

import pytest


@pytest.fixture
def shared_datadir(request, tmpdir):
    """Fixture to give access to a data directory.

    Taken from:
    https://github.com/gabrielcnr/pytest-datadir
    """

    original_shared_path = os.path.join(request.fspath.dirname, 'data')
    temp_path = pathlib.Path(str(tmpdir.join('data')))
    shutil.copytree(original_shared_path, str(temp_path))
    return temp_path
