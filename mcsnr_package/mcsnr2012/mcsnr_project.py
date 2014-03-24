from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from mcsnr2012.mcsnr_alchemy import Base

class MCSNRProject(object):

    def __init__(self, database_string, echo=False):
        self.metadata = Base.metadata
        self.engine = create_engine(database_string, echo=echo)
        self.metadata.bind = self.engine
        self.metadata.create_all()
        self.Session = sessionmaker(bind=self.engine)
        self.session = self.Session()
        self.conn = self.session.bind.connect()

