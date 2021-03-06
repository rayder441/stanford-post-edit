The Stanford Post-Editing Corpus
================================

Spence Green, Jeffrey Heer, and Christopher D. Manning

contact: http://nlp.stanford.edu // spence@spencegreen.com


This corpus contains a set of English source sentences translated to Arabic, French, and German by professional human translators. Each source sentence has 16 independent translations in each language. The translations were collected using a simple web interface with two translation conditions: unaided and post-edit. In the post-edit condition, suggested translations from Google Translate were provided. However, the translators were free to delete the suggested translation. The source sentences were presented in document order, but
subjects could not view the full document context. After submission of each translation, no further revision was allowed.

Subjects completed the experiment under time pressure. The interface displayed an idle timer that prohibited pauses longer than three minutes. The idle timer reset upon any keystroke in the target textbox. Upon expiration, it triggered submission of any entered text. The duration was chosen to allow reflection but to ensure completion of the experiment during a single session.

We recorded all keyboard, mouse, and browser events along with timestamps. The source tokens were also placed in separate <span> elements so that we could record hover events.


### SELECTION OF SOURCE MATERIALS

We chose English as the source language and Arabic, French, and German as the target languages. The target languages were selected based on canonical word order. Arabic is Verb-Subject-Object (VSO), French is SVO, and German is SOV. Verbs are salient linguistic elements that participate in many syntactic relations, so we wanted to control the position of this variable for cross-linguistic modeling. 

We selected four paragraphs from four English Wikipedia articles. We deemed two of the paragraphs "easy" in terms of lexical and syntactic features, and the other two "hard." Subjects saw one easy and one hard document in each translation condition. We selected the passages according to two criteria: "Featured" status---meaning that the text quality was rated highly by the Wikipedia community---in at least two languages and universality of the topic. For example, many American pop culture articles are featured, yet they often contain high concentrations of idiosyncratic lexical items. The four topics we selected were the 1896 Olympics (easy), the flag of Japan (easy), Schizophrenia (hard), and the infinite monkey theorem (hard).


### SELECTION OF SUBJECTS

For each target language, we hired 16 self-described "professional" translators on oDesk. Most were freelancers with at least a bachelor's degree. Three had Ph.Ds. We advertised the job at a fixed price of $0.085 per source token ($52 in total), a common rate for general text. However, we allowed translators to submit bids so that they felt fairly compensated. We did not negotiate, but the bids centered close to our target price: Arabic, $50.50; French, $52.32; German, $49.93. We also recorded user profile information included scores on relevant oDesk skills tests.


### CORPUS ORGANIZATION

    README.md         This file
    ar/               Arabic translations
    de/               German translations
    en-source/        English source documents
    fr/               French translations
    google-translate/ Google Translate output of the English source
    processed/        CSV files of the translation data
    subjects/         CSV files describing the human subjects

The translation directories contain the following files indexed by unique user id. For example, for Arabic translator #10:

    10.tgt.txt                The target segments
    10.tgt.meta.txt           Metadata for each target segment
    10.tgt.actionlog.txt      The Javascript action logs for each segment
    10.tgt.actionlog.txt.csv  Post-processed action logs

The columns "meta" in the meta file are:

    is_machine  "t" if the segment was generated by a machine
    date        Date when the translation was completed
    ui_id       1, un-aided condition; 2, post-edit
    is_valid    "t" if the idle timer expired
    src_len     # source tokens
    tgt_len     # target tokens

The action logs contain all Javascript events in the DOM Level 3 specification except abort, keyup, mousemove, and wheel. Events are delimited by "|" and the format is:

    event_name time css_id <payload>

Time is measured in milliseconds since the synthetic start event, which as time 0. Event-specific payloads include information like x/y mouse coordinates and the jQuery id of the keystroke that was pressed.


### FURTHER INFORMATION

Most subdirectories contain a README file that describes the directory contents.
