import pytest

def pytest_addoption(parser):
    parser.addoption("--username", action="store", help="NASA Earthdata username")
    parser.addoption("--password", action="store", help="NASA Earthdata password")
    parser.addoption("--webdav", action="store", help="PO.DAAC Drive WebDAV password")

@pytest.fixture
def username(request):
    """ Returns NASA Earthdata username """
    return request.config.getoption("--username")

@pytest.fixture
def password(request):
    """ Returns NASA Earthdata password """
    return request.config.getoption("--password")

@pytest.fixture
def webdav(request):
    """ Returns PO.DAAC Drive WebDAV password """
    return request.config.getoption("--webdav")
