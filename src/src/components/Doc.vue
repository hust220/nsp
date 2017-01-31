<template>
  <el-row class="doc">
    <el-col :span="5" class="sidebar">
      <div class="doc-nav" v-html="doc_nav"></div>
    </el-col>
    <el-col :span="19" :offset="5" class="main">
      <div class="doc-content" v-html="content"></div>
    </el-col>
  </el-row>
</template>

<script>
  import { bus } from '../bus.js'
  import marked from 'marked'

  export default {
    data() {
      return {
        doc_nav: '',
        content: 'hi'
      }
    },

    methods: {
      fetch_doc_nav() {
        var url = '/static/docs/nav.md'
        this.$http.get(url).then(response => {
          this.doc_nav = marked(response.body)
        }, response => {
          //
        })
      },

      fetch_content(item) {
        var url = '/static/docs/' + item + '.md'
        this.$http.get(url).then(response => {
          this.content = marked(response.body)
        }, response => {
          //
        })
      }
    },

    created() {
      this.fetch_doc_nav()
      var item = (this.$route.params.item ? this.$route.params.item : 'index')
      this.fetch_content(item)
    },

    watch: {
      '$route'(to, from) {
        var item = (to.params.item ? to.params.item : 'index')
        this.fetch_content(item)
      }
    }
  }
</script>

<style>
  .sidebar {
    position: fixed;
    font-size: 0.85em;
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
    margin: 5px;
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