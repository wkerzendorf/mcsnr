from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from mcsnr2012 import mcsnr_alchemy as malchemy
from mcsnr2012.mcsnr_alchemy_mixins import GeminiUtilDBMixin

class MCSNRProject(object):

    quick_search_classes = [malchemy.Candidate]

    quick_search_classes = [malchemy.Candidate]
    def __init__(self, database_string, echo=False):
        self.metadata = malchemy.Base.metadata
        self.engine = create_engine(database_string, echo=echo)
        self.metadata.bind = self.engine
        self.metadata.create_all()
        self.Session = sessionmaker(bind=self.engine)
        self.session = self.Session()
        self.conn = self.session.bind.connect()

    @property
    def quick_search_table_names(self):
        return [class_model.__tablename__ for class_model in self.quick_search_classes]

    @property
    def observed_snrs(self):
        return self.session.query(malchemy.SNRGeminiTarget).filter_by(observed=True)

    def set_geminiutil_session(self, geminiutil_session):
        GeminiUtilDBMixin.geminiutil_session = geminiutil_session
