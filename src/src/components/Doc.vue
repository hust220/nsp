<template>
  <el-row class="doc">
    <el-col :span="6" class="sidebar">
      <div class="doc-nav" v-html="doc_nav"></div>
    </el-col>
    <el-col :span="18" :offset="6" class="main">
      <div class="doc-content" v-html="content"></div>
    </el-col>
  </el-row>
</template>

<script>
  import { bus } from '../bus.js'
  import { ScrollToControl } from '../scroll.js'
  import docCache from '../data.js'
  import marked from 'marked'

  export default {
    data() {
      return {
        doc_nav: '',
        content: ''
      }
    },

    methods: {
      scroll(item) {
        if (item) {
          ScrollToControl(item)
        } else {
          window.scrollBy(0, -document.body.scrollTop)
        }
      },

      fetch_doc_nav() {
        var url = 'static/docs/nav.md'
        this.$http.get(url).then(response => {
          this.doc_nav = marked(response.body)
        }, response => {
          //
        })
      },

      fetch_doc_content() {
        var params = this.$route.params
        var theme = (params.theme ? params.theme : 'install')
        var item = (params.item ? params.item : '')
        if (docCache.theme && theme === docCache.theme) this.scroll(item)
        else if (theme in docCache) {
          this.content = docCache[theme]
          docCache.theme = theme
          this.scroll(item)
        } else {
          var url = 'static/docs/' + theme + '.md'
          this.$http.get(url).then(response => {
            this.content = marked(response.body)
            docCache[theme] = this.content
            this.scroll(item)
          }, response => {
            //
          })
          docCache.theme = theme
        }
      }
    },

    created() {
      this.fetch_doc_nav()
      this.fetch_doc_content()
    },

    watch: {
      '$route'(to, from) {
        this.fetch_doc_content()
      }
    },

    destroyed() {
      docCache.theme = ''
    }

  }
</script>

<style>
  .sidebar {
    position: fixed;
    font-size: 0.85em;
    overflow-y:scroll;
    height: 100%;
    padding-left: 50px;
    padding-bottom: 100px;
  }

  .sidebar ul {
    margin: 10px 15px;
    padding: 0px;
  }

  .sidebar li {
    list-style-type: none;
    margin: 10px 0px;
  }

  .doc-nav >ul {
    margin: 10px 0px;
  }

  .doc-content {
    margin: 35px;
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif, "Apple Color Emoji", "Segoe UI Emoji", "Segoe UI Symbol";
  }

  .doc-content pre {
    background-color: #f7f7f7;
    padding: 16px;
  }

  .doc-content code {
    font-family: Consolas, "Liberation Mono", Menlo, Courier, monospace;
  }

  .doc-content table {
    border-collapse: collapse;
  }

  .doc-content tr:nth-child(2n) {
    background-color: #f8f8f8;
  }

  .doc-content td {
    padding: 6px 13px;
    border: 1px solid #ddd;
  }

</style>